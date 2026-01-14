# Iterator of all filled chunks of a LongSequence, in order
struct LongSeqFullChunks
    chunks::Memory{UInt64}
    len::Int
end

Base.eltype(::Type{LongSeqFullChunks}) = UInt64
Base.length(itr::LongSeqFullChunks) = itr.len

function Base.iterate(itr::LongSeqFullChunks, state::Int = 1)
    # We simply emits all the chunks from the underlying memory directly
    state > length(itr) && return nothing
    return (@inbounds(itr.chunks[state]), state + 1)
end

# This iterator yields all the filled chunks of a ResizableSeq, or LongSubSeq, in order.
# Note that it does not yield them as stored in their Memory; instead, it yields them as if
# they were stored with an offset of zero.
# Hence, we a LongSeqFullChunks(a) and SubSeqChunks(b) emits identical elements if a == b 
struct SubSeqFullChunks
    chunks::Memory{UInt64}
    index_offset::Int # First chunk to read is index_offset + 1
    len::Int # length of iterator
    offset::Int # offset of first chunk to load
end

Base.eltype(::Type{SubSeqFullChunks}) = UInt64
Base.length(itr::SubSeqFullChunks) = itr.len

function Base.iterate(itr::SubSeqFullChunks, state::Int = 0)
    # Since the chunks may have an offset, we need to handle it.
    state ≥ itr.len && return nothing
    return if iszero(itr.offset)
        # If the offset is zrero, the layout is similar to a LongSequence,
        # and we can emit the chunks from the memory directly
        (@inbounds(itr.chunks[state + itr.index_offset + 1]), state + 1)
    else
        # Else, a whole chunk is composed of the upper parts of one chunk
        # and the lower parts of the next chunk.
        # Note that in this branch, we have already checked that:
        # * There is at least one full chunk to emit
        # * The offset is nonzero, so this chunk is split among to elements in memory
        # Therefore, we know we are inbounds
        left = @inbounds itr.chunks[state + itr.index_offset + 1]
        right = @inbounds itr.chunks[state + itr.index_offset + 2]
        chunk = right_shift(left, itr.offset) | left_shift(right, 64 - itr.offset)
        (chunk, state + 1)
    end
end

# This returns ((tail::UInt64, tail_bits::UInt8), itr::LongSeqFullChunks), where
# itr is an iterator of all filled chunks, in order,
# and `tail` represents the final non-filled chunk, with the lower `tail_bits` bit
# being coding.
function chunk_stream(seq::LongSequence)
    EMPTY = (zero(UInt64), 0x00)
    bps = bits_per_symbol(seq)
    iszero(bps) && return (EMPTY, LongSeqFullChunks(seq.chunks, 0))
    nbits = length(seq) * bps
    tail_bits = rem(nbits % UInt, UInt(64))
    return if iszero(tail_bits)
        (EMPTY, LongSeqFullChunks(seq.chunks, length(seq.chunks)))
    else
        mask = left_shift(UInt64(1), tail_bits) - UInt64(1)
        tail = (seq.chunks[end] & mask, tail_bits % UInt8)
        itr = LongSeqFullChunks(seq.chunks, length(seq.chunks) - 1)
        (tail, itr)
    end
end

# Same behaviour as the method with LongSequence, but yields SubSeqFullChunks
# as the iterator
function chunk_stream(seq::Union{ResizableSeq, LongSubSeq})
    EMPTY = (zero(UInt64), 0x00)
    bps = bits_per_symbol(seq)
    iszero(bps) && return (EMPTY, SubSeqFullChunks(seq.chunks, 0, 0, 0))
    nbits = length(seq) * bps
    (n_full_chunks, tail_bits) = divrem(nbits % UInt, UInt(64))
    n_full_chunks = n_full_chunks % Int
    tail_bits = tail_bits % Int
    symbols_per_chunk = div(64, bps)

    # Index offset: The zero-based offset where we load the first chunk from
    index_offset = div(seq.offset % UInt, symbols_per_chunk % UInt) % Int

    # This is the number of unused bottom bits in the first chunk we load
    chunk_offset = (((seq.offset * bps) % UInt) % UInt(64)) % Int

    # All these parameters are valid, no matter how the tail looks.
    itr = SubSeqFullChunks(seq.chunks, index_offset, n_full_chunks, chunk_offset)

    # Note that if the sequence is empty, this branch is hit
    iszero(tail_bits) && return (EMPTY, itr)

    # We now know the sequence is not empty, and the tail is not empty.
    # Cases:
    # ---###..          Tail is entirely within a chunk
    # -----### ##...... Tail is split along two chunks

    (last_chunk_idx_offset, last_offset) = divrem(((length(seq) + seq.offset - 1) * bps) % UInt, UInt(64))
    bits_in_last = (last_offset + bps) % Int
    tail = if bits_in_last ≥ tail_bits
        # Encompassed in one chunk
        chunk = @inbounds seq.chunks[(last_chunk_idx_offset + 1) % Int]
        right_shift(chunk, bits_in_last - tail_bits)
    else
        # Split along two chunks
        left = @inbounds seq.chunks[last_chunk_idx_offset % Int]
        right = @inbounds seq.chunks[(last_chunk_idx_offset + 1) % Int]
        bits_in_left_chunk = tail_bits - bits_in_last
        right_shift(left, 64 - bits_in_left_chunk) | left_shift(right, bits_in_left_chunk)
    end
    tail &= left_shift(UInt64(1), tail_bits) - UInt64(1)
    return ((tail, tail_bits % UInt8), itr)
end

# Returns (head, tail, body), where
# (a, b) = head::Tuple{UInt64, UInt8} is the first non-whole coding chunk with
#    only the lower b bits of a set.
# Tail is same as head, but the last non-whole coding chunk. Still only the lower
# bits are set
# Body is a SubArray of the memory containing only the whole chunks.
# If the head and the body are the same coding chunk, then the tail is zero.

# Here, if each symbol is 8 bits
#  ...##### ######## ######## ###.....
#     ----- ----------------- ---
#       |        |             | tail: (0x0000000000XXXXXX, 0x18)
#       |        | body: view(seq.chunks, 2:3)
#       | head: (0x000000XXXXXXXXXX, 0x28)
#
function ordered_parts(seq::LongSequence)
    EMPTY = (UInt64(0), 0x00)
    bps = bits_per_symbol(seq)
    if iszero(bps) || isempty(seq)
        return (EMPTY, EMPTY, view(seq.chunks, 1:0))
    end

    n_bits = (length(seq) * bps) % UInt
    (n_full_chunks, tail_bits) = divrem(n_bits, UInt(64))

    tail = if iszero(tail_bits)
        EMPTY
    else
        chunk = @inbounds seq.chunks[(n_full_chunks + 1) % Int]
        chunk &= left_shift(UInt64(1), tail_bits) - 1
        (chunk, tail_bits % UInt8)
    end

    # TODO: Could use memory views
    body = @inbounds view(seq.chunks, 1:(n_full_chunks % Int))

    return (EMPTY, tail, body)
end

ordered_parts(x::ResizableSeq) = ordered_parts(LongSubSeq(x))

function ordered_parts(x::LongSubSeq)
    EMPTY = (UInt64(0), 0x00)
    bps = bits_per_symbol(x)
    if iszero(bps) || isempty(x)
        return (EMPTY, EMPTY, view(x.chunks, 1:0))
    end

    (start_index, start_offset) = index_offset(x, 1)
    (past_index, past_offset) = index_offset(x, length(x) + 1)

    # Compute the body.
    # We skip the first chunk if the start offset is not zero
    # ..###### ########
    #   ^ Nonzero offset - skip first index

    #  ########
    #  ^ Zero offset - do not skip
    body_start = start_index + !iszero(start_offset)

    # If past_offset is zero (see below), then the final full index
    # is one before the index of past
    # ######## ........
    #          ^

    # If past_offset is not zero, then we have a tail, and so final index
    # is also one before index of past
    # ######## ###.....
    #             ^
    body_stop = past_index - 1

    # Get head.
    # ...#####
    #    ^ Start offset is not zero: We have a head

    # #######.
    #        ^ Start index and past index are the same: We also have a head
    head = if !iszero(start_offset) | (start_index == past_index)
        chunk = @inbounds x.chunks[start_index]
        # Mask out unused bits - n.b: If head and tail is the same, then
        # the length may be shorter than what the offset implies,
        # i.e. in the case of ...##...
        head_len = min(bps * length(x), 64 - start_offset) % UInt8
        chunk &= left_shift(UInt64(1), head_len) - 1
        (chunk, head_len)
    else
        EMPTY
    end

    # Get tail
    # ####....
    #     ^ Past offset is not zero: We have a head
    # ..##....
    # EXCEPT if start index and past index is the same, in which case the tail
    # is the head
    tail = if !iszero(past_offset) & (start_index != past_index)
        # Note: If past_offset is not zero, then past_index is the same index
        # as the last symbol
        chunk = @inbounds x.chunks[past_index]
        chunk &= left_shift(UInt64(1), past_offset) - 1
        (chunk, past_offset % UInt8)
    else
        EMPTY
    end

    # TODO: Could use memory views
    body = @inbounds view(x.chunks, body_start:body_stop)
    return (head, tail, body)
end

# This abstracts over the encoding of chunks from some ASCII data
# and into encoded chunks.
# This only emits full chunks. As such, the vector must have a length
# divisible by the number of symbols per chunk.
# This struct exists because
# * It enables reuse across various callsites that uses ASCIIEncoding
# * The validity check where we check for nothing in `try_ascii_encode`
#   is unrolled, and only triggers an error once per chunk.
#   this significantly improves encoding speed
struct ASCIIFullChunkIterator{A <: Alphabet, V <: AbstractVector{UInt8}}
    v::V

    # Nearly all this work happens at compile time, so let's inline it
    @inline function ASCIIFullChunkIterator{A, V}(v::V) where {A <: Alphabet, V <: AbstractVector{UInt8}}
        Base.require_one_based_indexing(v)
        bps = bits_per_symbol(A())
        iszero(bps) && error("Bug: Should never be called with zero-BPS alphabet")
        if !iszero((length(v) * bps) % 64)
            error("Vector length does not match one whole chunk")
        end
        return new{A, V}(v)
    end
end

function ASCIIFullChunkIterator(A::Alphabet, v::AbstractVector{UInt8})
    return ASCIIFullChunkIterator{typeof(A), typeof(v)}(v)
end

# Either the chunk, or the 1-based index of the first bad byte, and offending byte
Base.eltype(::Type{<:ASCIIFullChunkIterator}) = Union{UInt64, Tuple{Int, UInt8}}

# Since a bad byte can terminate it at any time
Base.length(::ASCIIFullChunkIterator{A}) where {A} = Base.SizeUnknown()

function Base.iterate(itr::ASCIIFullChunkIterator{A}, state::Int = 1) where {A}
    state > length(itr.v) && return nothing
    bps = bits_per_symbol(A())
    n_symbols = div(64, bps)
    bad = false
    chunk = zero(UInt64)
    shift = 0
    for i in 0:(n_symbols - 1)
        byte = @inbounds itr.v[state + i]
        enc = try_ascii_encode(A(), byte)
        if enc === nothing
            bad = true
            enc = zero(UInt64)
        end
        chunk |= left_shift(enc, shift)
        shift += bps
    end
    return if bad
        (get_bad_byte, typemax(Int))
    else
        (chunk, state + n_symbols)
    end
end

@noinline function get_bad_byte(itr::ASCIIFullChunkIterator{A}, state::Int) where {A}
    bps = bits_per_symbol(A())
    n_symbols = div(64, bps)
    for i in state:(state + n_symbols - 1)
        byte = @inbounds itr.v[state + i]
        enc = try_ascii_encode(A(), byte)
        enc === nothing && return (i, byte)
    end
    # The enclosing function is only called if we encounter a bad byte
    # when encoding an ASCII chunk, so we should reach this bad byte again
    # in the loop above, and exit
    error("Unreachable")
end

# Use ASCIIFullChunkIterator to construct full chunk, and this function to
# construct partial chunks. This returns the partial chunk, or else
# (index, bad_byte), if the byte is bad
function ascii_encode_partial_chunk_inbounds(
        A::Alphabet,
        v::AbstractVector{UInt8},
    )::Union{UInt64, Tuple{Int, UInt8}}
    bps = bits_per_symbol(A)

    # This function should only ever be called on non-whole chunks.
    # Since this check is cheap, let's do it.
    symbols_per_chunk = div(64, bps)
    length(v) < symbols_per_chunk || error("Unreachable")
    chunk = UInt64(0)
    shift = 0
    for i in 1:length(v)
        byte = @inbounds v[i]
        enc = try_ascii_encode(A, byte)
        enc === nothing && return (i, byte)
        chunk |= left_shift(enc, shift)
        shift += bps
    end
    return chunk
end
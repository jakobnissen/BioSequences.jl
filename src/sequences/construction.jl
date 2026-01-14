LongSequence(seq::BioSequence{A}) where {A} = LongSequence{A}(seq)

function LongSequence{A}(src) where {A}
    scheme = EncodingScheme(A(), typeof(src))
    src = to_decoding_source(scheme, src)
    return build_sequence(LongSequence{A}, scheme, src)
end

function build_sequence(::Type{LongSequence{A}}, ::EncodingScheme, src) where {A}
    len = UInt(length(src))::UInt
    seq = unsafe_undef_biosequence(LongSequence{A}, len)
    i = 0
    for element in src
        v = convert(eltype(LongSequence{A}), element)
        seq[(i += 1)] = v
    end
    (i % UInt) == len || throw(ArgumentError("Too few elements in source"))
    return seq
end

function build_sequence(::Type{LongSequence{A}}, ::FourToTwo, src::HeapSeq) where A
    # It is more efficient to scan in one pass, then build the sequence in another,
    # because this allows each of the loops to be far more efficient, and possibly SIMD.
    bad_idx = findfirst(!iscertain, src)
    isnothing(bad_idx) || throw(EncodeError(A(), src[bad_idx], bad_idx))
    seq = unsafe_undef_biosequence(LongSequence{A}, length(src) % UInt)
    ((tail, _), itr) = chunk_stream(src)



    # TODO
    # 
end

# TODO: Four to two
# TODO: Two to four

function build_sequence(::Type{LongSequence{A}}, ::CopyableEncoding, seq::LongSequence) where {A}
    return unsafe_new_seq(LongSequence{A}, copy(seq.chunks), length(seq))
end

function build_sequence(::Type{LongSequence{A}}, ::CopyableEncoding, seq::Union{ResizableSeq, LongSubSeq}) where {A}
    ((head, nhead), itr) = iter_chunks(seq)
    mem = Memory{UInt64}(undef, length(itr) + !iszero(nhead))
    if !iszero(nhead)
        mem[end] = head
    end
    for (i, chunk) in enumerate(itr)
        mem[i] = chunk
    end
    return unsafe_new_seq(LongSequence{A}, mem, length(seq))
end

function build_sequence(::Type{LongSequence{A}}, ::ASCIIEncoding, src::AbstractVector{UInt8}) where {A}
    len = UInt(length(src))::UInt
    seq = unsafe_undef_biosequence(LongSequence{A}, len)
    bps = bits_per_symbol(A())
    iszero(bps) && return seq

    # Construct full chunks
    symbols_per_chunk = div(64, bps)
    (body_chunks, tail_symbols) = divrem(length(src) % UInt, symbols_per_chunk)
    body_chunks = body_chunks % Int
    if body_chunks > 0
        vw = @inbounds view(src, 1:(body_chunks * symbols_per_chunk))
        for (i, chunk) in enumerate(ASCIIFullChunkIterator(A(), vw))
            if chunk isa UInt64
                seq.chunks[i] = chunk
            else
                (index, byte) = chunk
                throw(EncodeError(A(), byte, index))
            end
        end
    end

    # Construct tail
    if tail_symbols > 0
        vw = @inbounds view(src, (body_chunks * symbols_per_chunk + 1):length(src))
        chunk = ascii_encode_partial_chunk_inbounds(A(), vw)
        if chunk isa UInt64
            seq.chunks[end] = chunk
        else
            (index, byte) = chunk
            throw(EncodeError(A(), byte, index))
        end
    end

    return seq
end

# LongSubSeq constructors
# Note: These are views, and only supports the three BioSequence types defined in this package,
# which has a memory layout compatible with LongSubSeq.
# We don't do any copying constructors, since this is a view.

LongSubSeq(src::Union{LongSequence{A}, ResizableSeq{A}, LongSubSeq{A}}) where {A} = LongSubSeq{A}(src)

function LongSubSeq{A}(src) where {A}
    scheme = EncodingScheme(A(), typeof(src))
    src = to_decoding_source(scheme, src)
    return build_sequence(LongSubSeq{A}, scheme, src)
end

function build_sequence(::Type{LongSubSeq{A}}, ::CopyableEncoding, src::Union{LongSubSeq, ResizableSeq}) where {A}
    return unsafe_new_seq(LongSubSeq{A}, src.chunks, src.offset, src.len)
end

function build_sequence(::Type{LongSubSeq{A}}, ::CopyableEncoding, src::LongSequence{A}) where {A}
    return unsafe_new_seq(LongSubSeq{A}, src.chunks, 0, src.len)
end

function Base.view(src::LongSequence{A}, idx::AbstractUnitRange) where {A}
    isempty(idx) && return unsafe_new_seq(LongSubSeq{A}, src.chunks, 0, 0)
    @boundscheck checkbounds(src, idx)
    return unsafe_new_seq(LongSubSeq{A}, src.chunks, Int(first(idx) - 1)::Int, Int(length(idx))::Int)
end

function Base.view(src::Union{LongSubSeq{A}, ResizableSeq{A}}, idx::AbstractUnitRange) where {A}
    isempty(idx) && return unsafe_new_seq(LongSubSeq{A}, src.chunks, 0, 0)
    @boundscheck checkbounds(src, idx)
    offset = src.offset + Int(first(idx))::Int - 1
    return unsafe_new_seq(LongSubSeq{A}, src.chunks, offset, Int(length(idx))::Int)
end

# TODO: Unsafe interpretation of memory as LongSubSeq

# TODO: Resizable, all the same as LongSequence
ResizableSeq(seq::BioSequence{A}) where {A} = ResizableSeq{A}(seq)

function ResizableSeq{A}(src) where {A}
    scheme = EncodingScheme(A(), typeof(src))
    src = to_decoding_source(scheme, src)
    return build_sequence(ResizableSeq{A}, scheme, src)
end

# When building a new sequence, offset is guaranteed to be zero, and when the offset is zero,
# the memory layout is identical to LongSequence, so we can just use the LongSequence
# constructors internally
function build_sequence(::Type{ResizableSeq{A}}, E::EncodingScheme, src) where {A}
    long_seq = build_sequence(LongSequence{A}, E, src)
    unsafe_new_seq(ResizableSeq{A}, long_seq.chunks, 0, length(long_seq))
end
# Same as `get_parts`, except does not return a head, because
# LongSequence never has a head
function get_longseq_parts(seq::LongSequence)
    bps = bits_per_symbol(seq)
    if iszero(bps)
        return ((UInt64(0), 0x00), view(seq.chunks, 1:0))
    end

    n_bits = (length(seq) * bps) % UInt
    (n_full_chunks, tail_bits) = divrem(n_bits, UInt(64))

    tail = if iszero(tail_bits)
        (UInt64(0), 0x00)
    else
        chunk = @inbounds seq.chunks[(n_full_chunks + 1) % Int]
        (chunk, tail_bits % UInt8)
    end

    # TODO: Could use memory views
    body = @inbounds view(seq.chunks, 1:(n_full_chunks % Int))
    return (tail, body)
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
function get_parts(seq::LongSequence)
    (head, body) = get_longseq_parts(seq)
    tail = (UInt64(0), 0x00)
    return (head, tail, body)
end

get_parts(x::ResizableSeq) = get_parts(LongSubSeq(x))

function get_parts(x::LongSubSeq)
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

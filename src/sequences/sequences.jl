function unsafe_new_seq end

"""
    LongSequence{A <: Alphabet} <: BioSequence{A}

General-purpose `BioSequence`. This type is mutable and variable-length, but
not resizable. It is be preferred for most use cases.

See also: [`BioSequence`](@ref), [`Alphabet`](@ref), [`LongSubSeq`](@ref), [`ResizableSeq`](@ref)

# Extended help
The following table summarizes common LongSequence types that have been given
aliases for convenience.

| Type                                | Symbol type | Type alias   |
| :---------------------------------- | :---------- | :----------- |
| `LongSequence{DNAAlphabet{N}}`      | `DNA`       | `LongDNA{N}` |
| `LongSequence{RNAAlphabet{N}}`      | `RNA`       | `LongRNA{N}` |
| `LongSequence{AminoAcidAlphabet}`   | `AminoAcid` | `LongAA`     |

The `LongDNA` and `LongRNA` aliases use a four bit alphabet (e.g. `DNAAlphabet{4}`).
"""
struct LongSequence{A <: Alphabet} <: BioSequence{A}
    # Sequences chunks are stored from first to last index, and
    # within each chunk, from lowest to highest bit.
    # The length is always cld(64, bps), where bps is the
    # bits per symbol of A, and chunks is always empty
    # if bps is zero.
    # Unused bits in the last element may have any value,
    # and may be freely mutated. When read, they should be masked out.
    chunks::Memory{UInt64}

    # Number of symbols. Non-negative.
    len::Int

    global function unsafe_new_seq(
            ::Type{LongSequence{A}},
            chunks::Memory{UInt64},
            len::Int
        ) where {A <: Alphabet}
        return new{A}(chunks, len)
    end
end

const LongDNA{N} = LongSequence{DNAAlphabet{N}}
const LongRNA{N} = LongSequence{RNAAlphabet{N}}
const LongAA = LongSequence{AminoAcidAlphabet}

"""
    ResizableSeq{A <: Alphabet} <: BioSequence{A}

Mutable, resizable `BioSequence` type. This type is used over `LongSequence`
for its ability to be resized.

See also: [`BioSequence`](@ref), [`Alphabet`](@ref), [`LongSequence`](@ref), [`LongSubSeq`](@ref)
"""
mutable struct ResizableSeq{A <: Alphabet} <: BioSequence{A}
    # Same encoding as LongSequence, but length may be longer
    # than the bare minimum necessity
    chunks::Memory{UInt64}

    # The first offset * bps bits are unused, and follow the semantics
    # of unused bits in LongSequence
    offset::Int

    # Number of symbols. Non-negative.
    len::Int

    global function unsafe_new_seq(
            ::Type{ResizableSeq{A}},
            chunks::Memory{UInt64},
            offset::Int,
            len::Int
        ) where {A <: Alphabet}
        return new{A}(chunks, offset, len)
    end
end

"""
    LongSubSeq{A <: Alphabet} <: BioSequence{A}

This sequence is similar to `LongSequence`, but is used as a view into an existing
`LongSequence`, `ResizableSeq`.
It can also be used to interpret a `Memory{UInt64}` as a `BioSequence`.

See also: [`BioSequence`](@ref), [`Alphabet`](@ref), [`LongSequence`](@ref), [`ResizableSeq`](@ref)
"""
struct LongSubSeq{A <: Alphabet} <: BioSequence{A}
    # Fields mean the exact same as `ResizableSeq`.
    chunks::Memory{UInt64}
    offset::Int
    len::Int

    global function unsafe_new_seq(
            ::Type{LongSubSeq{A}},
            chunks::Memory{UInt64},
            offset::Int,
            len::Int
        ) where {A <: Alphabet}
        return new{A}(chunks, offset, len)
    end
end

@inline function left_shift(x::Unsigned, n::Integer)
    return x << (n & ((sizeof(x) * 8) - 1))
end

@inline function right_shift(x::Unsigned, n::Integer)
    return x >>> (n & ((sizeof(x) * 8) - 1))
end

get_bitmask(::BitsPerSymbol{N}) where {N} = UInt64(1) << N - UInt64(1)
get_bitmask(::BioSequence{A}) where {A} = get_bitmask(BitsPerSymbol(A()))

function get_min_chunks(T::Type{<:Union{LongSequence, ResizableSeq, LongSubSeq}}, len::UInt)
    bps = bits_per_symbol(T)
    return iszero(bps) ? 0 : cld(len * (bps % UInt), UInt64(64)) % Int
end

function index_offset(seq::LongSequence, i::Int)
    bps = bits_per_symbol(seq)
    # With zero BPS, we should never index any memory, so index zero
    # is indeed correct.
    iszero(bps) && return (0, 0)
    bitoffset = ((i - 1) % UInt) * (bps % UInt)
    (i, o) = divrem(bitoffset, UInt(64))
    return ((i % Int) + 1, o % Int)
end

function index_offset(seq::Union{ResizableSeq, LongSubSeq}, i::Int)
    bps = bits_per_symbol(seq)
    iszero(bps) && return (0, 0)
    bitoffset = ((i + seq.offset - 1) % UInt) * (bps % UInt)
    (i, o) = divrem(bitoffset, UInt(64))
    return ((i % Int) + 1, o % Int)
end

Base.length(x::Union{LongSequence, LongSubSeq, ResizableSeq}) = x.len

function Base.copy(x::LongSequence)
    return unsafe_new_seq(typeof(x), copy(x.chunks), x.len)
end

function Base.copy(x::Union{ResizableSeq, LongSubSeq})
    mem = Memory{UInt8}(undef, get_min_chunks(typeof(x), x.len))
    seq = unsafe_new_seq(typeof(x), mem, 0, x.len % Int)
    return @inbounds copy!(seq, x)
end

function unsafe_get_encoding(x::LongSequence, i::Int)
    iszero(bits_per_symbol(x)) && return zero(UInt64)
    (i, o) = index_offset(x, i)
    chunk = @inbounds x.chunks[i]
    return right_shift(chunk, o) & get_bitmask(x)
end

function unsafe_set_encoding!(x::LongSequence, enc::UInt64, i::Int)
    iszero(bits_per_symbol(x)) && return x
    (i, o) = index_offset(x, i)
    chunk = @inbounds x.chunks[i]
    chunk &= ~left_shift(get_bitmask(x), o)
    chunk |= left_shift(enc, o)
    @inbounds x.chunks[i] = chunk
    return x
end

function unsafe_get_encoding(x::Union{ResizableSeq, LongSubSeq}, i::Int)
    iszero(bits_per_symbol(x)) && return zero(UInt64)
    (i, o) = index_offset(x, i)
    chunk = @inbounds x.chunks[i]
    return right_shift(chunk, o) & get_bitmask(x)
end

function unsafe_set_encoding!(x::Union{ResizableSeq, LongSubSeq}, enc::UInt64, i::Int)
    iszero(bits_per_symbol(x)) && return x
    (i, o) = index_offset(x, i)
    chunk = @inbounds x.chunks[i]
    chunk &= ~left_shift(get_bitmask(x), o)
    chunk |= left_shift(enc, o)
    @inbounds x.chunks[i] = chunk
    return x
end

function unsafe_undef_biosequence(T::Type{<:LongSequence}, len::UInt64)
    return unsafe_new_seq(T, Memory{UInt64}(undef, get_min_chunks(T, len)), len % Int)
end

function unsafe_undef_biosequence(T::Type{<:Union{ResizableSeq, LongSubSeq}}, len::UInt64)
    return unsafe_new_seq(T, Memory{UInt64}(undef, get_min_chunks(T, len)), 0, len % Int)
end

function truncate!(x::ResizableSeq, len::UInt)
    @boundscheck checkbounds(x, len % Int)
    return unsafe_resize!(x, len)
end

function unsafe_resize!(x::ResizableSeq, len::UInt)
    x.len = len % Int
    return x
end

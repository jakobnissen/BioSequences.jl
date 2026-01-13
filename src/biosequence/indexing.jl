###
### Indexing
###
###
### Indexing methods for mutable biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.iterate(seq::BioSequence, i::Int = 1)
    (i % UInt) - 1 < (lastindex(seq) % UInt) ? (@inbounds seq[i], i + 1) : nothing
end

## Bounds checking
function Base.checkbounds(x::BioSequence, i::Integer)
    firstindex(x) ≤ i ≤ lastindex(x) || throw(BoundsError(x, i))
end

function Base.checkbounds(x::BioSequence, locs::AbstractVector{Bool})
    length(x) == length(locs) || throw(BoundsError(x, lastindex(locs)))
end

@inline function Base.checkbounds(seq::BioSequence, locs::AbstractVector)
    for i in locs
        checkbounds(seq, i)
    end
    return true
end

@inline function Base.checkbounds(seq::BioSequence, range::UnitRange)
    if !isempty(range) && (first(range) < 1 || last(range) > length(seq))
        throw(BoundsError(seq, range))
    end
end

## Getindex
function Base.getindex(x::BioSequence, i::Integer)
    i = Int(i)::Int
    @boundscheck checkbounds(x, i)
    data = unsafe_get_encoding(x, i)
    return decode(Alphabet(x), data)
end

function inbounds_copy_element!(dst::BioSequence, di::Int, src, si::Int)
    inbounds_copy_element!(EncodingScheme(Alphabet(dst), typeof(src)), di, src, si)
end

function inbounds_copy_element!(
        ::EncodingScheme,
        dst::BioSequence,
        di::Int,
        src,
        si::Int
    )
    symbol = @inbounds src[si]
    symbolT = convert(eltype(dst), symbol)
    @inbounds dst[di] = symbolT
end

function inbounds_copy_element!(
        ::CopyableEncoding,
        dst::BioSequence,
        di::Int,
        src::BioSequence,
        si::Int
    )
    enc = unsafe_get_encoding(src, si)
    unsafe_set_encoding!(dst, enc, di)
end

function inbounds_copy_element!(
        ::ASCIIEncoding,
        dst::BioSequence,
        di::Int,
        src,
        si::Int
    )
    byte = src[si]::UInt8
    enc = try_ascii_encode(Alphabet(dst), byte)
    isnothing(enc) && throw_byte_encoding(Alphabet(dst), byte)
    unsafe_set_encoding!(dst, enc, di)
end

@noinline function throw_byte_encoding(A::Alphabet, byte::UInt8)
    throw(EncodeError(A, Char(byte)))
end

function inbounds_copy_element!(
        ::TwoToFour,
        dst::BioSequence,
        di::Int,
        src::BioSequence,
        si::Int
    )
    enc = unsafe_get_encoding(src, si)
    enc = UInt64(1) << enc
    unsafe_set_encoding!(dst, enc, di)
end

function inbounds_copy_element!(
        ::FourToTwo,
        dst::BioSequence,
        di::Int,
        src::BioSequence,
        si::Int
    )
    enc = unsafe_get_encoding(src, si)
    if count_ones(enc % UInt8) != 1
        throw_fourtotwo_encoding(Alphabet(dst), Alphabet(src), enc)
    end
    enc = trailing_zeros(enc) % UInt64
    unsafe_set_encoding!(dst, enc, di)
end

@noinline function throw_fourtotwo_encoding(dstA::Alphabet, srcA::Alphabet, enc::UInt64)
    symbol = decode(srcA, enc)
    throw(EncodeError(dstA, symbol))
end

# Generic method: Create undefined sequence, then set each element
function Base.getindex(x::BioSequence, bools::AbstractVector{Bool})
    @boundscheck checkbounds(x, bools)
    res = unsafe_undef_biosequence(typeof(x), count(bools) % UInt)
    di = 0
    for (si, bool) in enumerate(bools)
        if bool
            di += 1
            inbounds_copy_element!(res, di, x, si)
        end
    end
    # For safety. We don't need to, but this check is cheap.
    di == length(res) || error("Element count does not match")
    return res
end

function Base.getindex(x::BioSequence, idx::AbstractVector{<:Integer})
    @boundscheck checkbounds(x, idx)
    res = unsafe_undef_biosequence(typeof(x), length(idx) % UInt)
    for (di, si) in enumerate(idx)
        inbounds_copy_element!(res, di, x, si)
    end
    return res
end

Base.getindex(x::BioSequence, ::Colon) = copy(x)

function Base.getindex(x::BioSequence, idx::AbstractUnitRange)
    res = unsafe_undef_biosequence(typeof(x), length(idx) % UInt)
    vw = view(x, idx)
    copy!(res, vw)
    return res
end

## Setindex
function Base.setindex!(x::BioSequence, v, i::Integer)
    i = Int(i)::Int
    @boundscheck checkbounds(x, i)
    vT = convert(eltype(typeof(x)), v)
    enc = encode(Alphabet(x), vT)
    unsafe_set_encoding!(x, enc, i)
end

function Base.setindex!(seq::BioSequence, x, locs::AbstractVector{<:Integer})
    @boundscheck checkbounds(seq, locs)
    @boundscheck if length(x) != length(locs)
        throw(DimensionMismatch("Attempt to assign $(length(x)) values to $(length(locs)) destinations"))
    end
    for (si, di) in enumerate(locs)
        inbounds_copy_element!(seq, Int(di)::Int, x, si)
    end
    return seq
end

function Base.setindex!(seq::BioSequence, x, locs::AbstractVector{Bool})
    @boundscheck checkbounds(seq, locs)
    n = count(locs)
    @boundscheck if length(x) != n
        throw(DimensionMismatch("Attempt to assign $(length(x)) values to $n destinations"))
    end
    si = 0
    for (di, bool) in enumerate(locs)
        if bool
            si += 1
            inbounds_copy_element!(seq, di, x, si)
        end
    end
    return seq
end

function Base.setindex!(seq::BioSequence, x, ::Colon)
    return setindex!(seq, x, 1:lastindex(seq))
end

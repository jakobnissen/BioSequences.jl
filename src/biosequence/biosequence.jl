###
### An abstract biological sequence type.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    abstract type BioSequence{A <: Alphabet}

A `BioSequence` stores `BioSymbol`s in a linear container, to represent polymer
molecules such as DNA, RNA, or proteins.
The default `BioSequence` type is `LongSequence{A}`.

See also: [`Alphabet`](@ref), [`LongSequence`](@ref)

# Extended help
The associated `Alphabet` type controls which `BioSymbol`s can be stored in the
`BioSequence`, and how these symbols are encoded to a binary encoding.
The `BioSequence` instance then handles the storage and retrieval of these
binary encodings.

!!! warning
    Decoding may assume the encoding is invalid, so a `BioSequence` may
    **never** give an invalid encoding to a user.
    This implies that `BioSequence`s that allow uninitialized data must
    either do so through `unsafe_` functions, or validate encodings before
    passing them onto decoding functions.

Subtypes of `BioSequence` must follow the below semantics:
* They should be indexable by `Int`, with indices `Base.OneTo(length(x))`
* They have an `Int` length, and iterating them is equivalent to indexing
  each of their elements in order, similar to a `Vector`

Subtypes `S <: BioSequence` must implement:
* `Base.length(::S)::Int`
* `unsafe_get_encoding(::S, ::Int)::UInt64`
* Instances of `S` must be constructable from any iterable with `length` defined and
  with a known, compatible element type.
* `Base.copy(::S)::S`

Mutable `BioSequence`s ought to define
* `unsafe_set_encoding!(::S, ::UInt64, ::Int)`
* `truncate!`
* `unsafe_resize!`
* `view`, or `fill_resize!`
* `fill_encoding!`
"""
abstract type BioSequence{A <: Alphabet} end

# Specific biosequences
"An alias for `BioSequence{<:NucleicAcidAlphabet}`"
const NucSeq = BioSequence{<:NucleicAcidAlphabet}

"An alias for `BioSequence{DNAAlphabet{N}}`"
const DNASeq{N} = BioSequence{DNAAlphabet{N}}

"An alias for `BioSequence{RNAAlphabet{N}}`"
const RNASeq{N} = BioSequence{RNAAlphabet{N}}

"An alias for `BioSequence{AminoAcidAlphabet}`"
const AASeq = BioSequence{AminoAcidAlphabet}

Base.eachindex(x::BioSequence) = Base.OneTo(length(x))
Base.firstindex(::BioSequence) = 1
Base.lastindex(x::BioSequence) = length(x)
Base.keys(seq::BioSequence) = eachindex(seq)
Base.nextind(::BioSequence, i::Integer) = Int(i) + 1
Base.prevind(::BioSequence, i::Integer) = Int(i) - 1
Base.eltype(::Type{<:BioSequence{A}}) where {A <: Alphabet} = eltype(A)

Alphabet(::Type{<:BioSequence{A}}) where {A <: Alphabet} = A()
Alphabet(x::BioSequence) = Alphabet(typeof(x))

Base.isempty(x::BioSequence) = iszero(length(x))
Base.empty(::Type{T}) where {T <: BioSequence} = T(eltype(T)[])
Base.empty(x::BioSequence) = empty(typeof(x))

BitsPerSymbol(x::BioSequence) = BitsPerSymbol(Alphabet(typeof(x)))
bits_per_symbol(::Type{T}) where {T <: BioSequence} = bits_per_symbol(Alphabet(T))
bits_per_symbol(x::BioSequence) = bits_per_symbol(typeof(x))

"""
    unsafe_get_encoding(s::BioSequence, i::Int)::UInt64

Extract the encoding of `s` at index `i`.
Subtypes of `BioSequence` must ensure this function never returns
an invalid encoding for the alphabet of `s`.

!!! warning
    This function does not boundscheck, and may result in undefined
    behaviour if `i` is out of bounds.

See also: [`unsafe_set_encoding!`]
"""
function unsafe_get_encoding end

"""
    unsafe_set_encoding!(s::BioSequence, enc::UInt64, i::Int) -> s

Store the encoding `enc` at index `i` in `s`. This operation is what
is responsible for `setindex!` of biosequences.

!!! warning
    This function must be called with an inbounds `i`.
    Furthermore, it assumes `enc` is a valid encoding in the alphabet
    of `s`. If either condition is violated, undefined behaviour may result.

See also: [`unsafe_get_encoding`]
"""
function unsafe_set_encoding! end

"""
    unsafe_undef_biosequence(T::Type{<:BioSequence}, len::UInt)::T

Create an unitialized `BioSequence` of type `T` and length `len`.

!!! warning
  This function is unsafe, because `BioSequence`s must never load
  symbols from invalid encodings.
  Hence, valid use of this function requires the caller to fill in
  all used coding bits of the sequence with a valid value.
"""
function unsafe_undef_biosequence end

"""
    unsafe_resize!(x::BioSequence, len::UInt) -> x

Resize `x` to length `len`. If `len` is shorter than the current length,
truncate the sequence. If `len` is larger, add new elements, which may be
uninitialized.

!!! warning
  This function is unsafe, because `BioSequence`s must never load
  symbols from invalid encodings.
  Hence, valid use of this function requires the caller to fill in
  all used coding bits of the sequence with a valid value.
"""
function unsafe_resize! end

"""
    fill_resize!(x::BioSequence, sym::BioSymbol, len::UInt)

Resize `x` to length `len`. If `len` is shorter than the current length,
truncate the sequence. If `len` is larger, push `sym` elements to the end
until length matches.
"""
function fill_resize!(x::BioSequence, sym::BioSymbol, len::UInt)
    oldlen = (length(x) % UInt)
    len == oldlen && return x
    len < oldlen && return @inbounds truncate!(x, len)
    enc = encode(Alphabet(x), sym)::UInt64
    unsafe_resize!(x, len)
    vw = @inbounds view(x, ((oldlen + 1) % Int):(len % Int))
    fill_encoding!(vw, enc)
    return x
end


function fill!(x::BioSequence, sym::BioSymbol)
    enc = encode(Alphabet(x), sym)::UInt64
    fill_encoding!(x, enc)
end

function fill_encoding!(x::BioSequence, enc::UInt64)
    for i in eachindex(x)
        unsafe_set_encoding!(x, enc, i)
    end
    return x
end

"""
    truncate!(x::BioSequence, len::UInt) -> x

Truncate `x` to length `len`. If `len` is larger than the current length,
throw a `BoundsError`
"""
function truncate! end
#=



function (::Type{S})(seq::BioSequence) where {S <: AbstractString}
    _string(S, seq, codetype(Alphabet(seq)))
end

Base.LazyString(seq::BioSequence) = LazyString(string(seq))

function _string(::Type{S}, seq::BioSequence, ::AlphabetCode) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end

function _string(::Type{String}, seq::BioSequence, ::AsciiAlphabet)
    String([stringbyte(s) for s in seq])
end

"""
    has_interface(::Type{BioSequence}, ::T, syms::Vector, mutable::Bool, compat::Bool=true)

Check if type `T` conforms to the `BioSequence` interface. A `T` is constructed from the vector
of element types `syms` which must not be empty.
If the `mutable` flag is set, also check the mutable interface.
If the `compat` flag is set, check for compatibility with existing alphabets.
"""
function has_interface(
    ::Type{BioSequence},
    ::Type{T},
    syms::Vector,
    mutable::Bool,
    compat::Bool=true
) where {T <: BioSequence}
    try
        isempty(syms) && error("Vector syms must not be empty")
        first(syms) isa eltype(T) || error("Vector is of wrong element type")
        seq = T((i for i in syms))
        length(seq) == length(syms) || return false
        eachindex(seq) === Base.OneTo(length(seq)) || return false
        E = encoded_data_eltype(T)
        e = extract_encoded_element(seq, 1)
        e isa E || return false
        (!compat || E == UInt) || return false
        copy(seq) isa typeof(seq) || return false
        if mutable
            encoded_setindex!(seq, e, 1)
            T(undef, 5) isa T || return false
            isempty(resize!(seq, 0)) || return false
        end
    catch error
        error isa MethodError && return false
        rethrow(error)
    end
    return true
end

function Base.similar(seq::BioSequence, len::Integer=length(seq))
    return typeof(seq)(undef, len)
end

# Fast path for iterables we know are stateless
function join!(seq::BioSequence, it::Union{Vector, Tuple, Set})
    _join!(resize!(seq, reduce((a, b) -> a + joinlen(b), it, init=0)), it, Val(true))
end

"""
    join!(seq::BioSequence, iter)

Concatenate all biosequences/biosymbols in `iter` into `seq`, resizing it to fit.

# Examples
```
julia> join(LongDNA(), [dna"TAG", dna"AAC"])
6nt DNA Sequence:
TAGAAC
```

see also [`join`](@ref)
"""
join!(seq::BioSequence, it) = _join!(seq, it, Val(false))

# B is whether the size of the destination seq is already
# known to be the final size
function _join!(seq::BioSequence, it, ::Val{B}) where B
    len = 0
    oldlen = length(seq)
    for i in it
        pluslen = joinlen(i)
        if !B && oldlen < (len + pluslen)
            resize!(seq, len + pluslen)
        end
        if i isa BioSymbol
            seq[len + 1] = i
        else
            copyto!(seq, len + 1, i, 1, length(i))
        end
        len += pluslen
    end
    seq
end

"""
    join(::Type{T <: BioSequence}, seqs)

Concatenate all the biosequences/biosymbols in `seqs` to a biosequence of type `T`.

# Examples
```
julia> join(LongDNA, [dna"TAG", dna"AAC"])
6nt DNA Sequence:
TAGAAC
```

see also [`join!`](@ref)
"""
function Base.join(::Type{T}, it::Union{Vector, Tuple, Set}) where {T <: BioSequence}
    _join!(T(undef, reduce((a, b) -> a + joinlen(b), it, init=0)), it, Val(true))
end

# length is intentionally not implemented for BioSymbol
joinlen(x::Union{BioSequence, BioSymbol}) = x isa BioSymbol ? 1 : length(x)

function Base.join(::Type{T}, it) where {T <: BioSequence}
    _join!(empty(T), it, Val(false))
end

Base.repeat(chunk::BioSequence, n::Integer) = join(typeof(chunk), (chunk for i in 1:n))
Base.:^(x::BioSequence, n::Integer) = repeat(x, n)

# Concatenation and Base.repeat operators
function Base.:*(fst::BioSequence, rest::BioSequence...)
    T = typeof(fst)
    join(T, (fst, rest...))
end


# The generic functions for any BioSequence...
include("indexing.jl")
include("predicates.jl")
include("find.jl")
include("printing.jl")
include("transformations.jl")
include("copying.jl")

=#

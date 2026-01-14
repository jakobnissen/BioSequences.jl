###
### Alphabet
###
###
### Alphabets of biological symbols.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    abstract type Alphabet

`Alphabet` is the most important type trait for `BioSequence`. An `Alphabet`
represents a set of biological symbols that can be encoded by a sequence,
e.g. A, C, G and T for a DNA Alphabet that requires only 2 bits
to represent each symbol.

`Alphabet` subtypes control the encoding compatible `BioSymbol`s into a
binary `UInt64` encoding, and the inverse decoding.
Subtypes of `BioSequence` then control how to store and retrieve the
encodings.

See also: [`BioSequence`](@ref)

# Extended help
* Subtypes `A <: Alphabet` are singleton structs constructable by `A()`.
* Alphabets define a *finite* set of biological symbols.
* The alphabet controls the encoding from a `BioSymbol` of the alphabet's element type,
  to an `UInt64`, as well as the decoding, the inverse process.
* Encoding is bijective: A symbol may only encode to one encoding, and one encoding
  only decodes to one symbol
* An `Alphabet`'s `encode` method must not produce invalid data.

### Required methods
Every subtype `A` of `Alphabet` must implement:
* `Base.eltype(::Type{A})::Type{S <: BioSymbol}`.
* `Base.iterate(::A, state)` iterate over all instances of `S` encodable in `A`
* `Base.length(::A)` the number of encodable elements in `A`
* `BioSequences.tryencode(::A, ::S)::UInt64 where {S <: BioSymbol}` encodes a symbol
  to the binary encoding.
* `BioSequences.decode(::A, ::UInt64)::S` decodes a binary encoding back to a symbol
* `BitsPerSymbol(::A)::BitsPerSymbol{N}`, where the `N` must be zero
  or a power of two in [1, 2, 4, 8, 16, 32, 64].

### Optional methods
* `AlphabetCode(::A)` for increased printing/writing efficiency
* `try_ascii_encode(::A, ::UInt8)`
* `ascii_decode(::A, ::UInt64)`
* `iscomplete(::A)`
"""
abstract type Alphabet end

"""
    struct BitsPerSymbol{N}

A trait object specifying the number of bits it takes to encode a `BioSymbol` in an `Alphabet`.
Alphabets `A` should implement `BitsPerSymbol(::A)::BitsPerSymbol`.
For compatibility with existing BioSequences, the number of bits should be a power of two
between 0 and 64, both inclusive.

See also: [`Alphabet`](@ref)

# Examples
```jldoctest
julia> BioSequences.BitsPerSymbol{8}()
BioSequences.BitsPerSymbol{8}()

julia> BioSequences.BitsPerSymbol(AminoAcidAlphabet())
BioSequences.BitsPerSymbol{8}()

julia> BioSequences.BitsPerSymbol(RNAAlphabet{2}())
BioSequences.BitsPerSymbol{2}()
```
"""
struct BitsPerSymbol{N} end

# Get the number of bits per symbol as an integer
bits_per_symbol(A::Alphabet) = bits_per_symbol(BitsPerSymbol(A))

bits_per_symbol(::BitsPerSymbol{N}) where {N} = N

"""
    iscomplete(::Alphabet)::Bool

Implement `iscomplete(::A) = true` for new subtypes `A <: Alphabet`
if `A` encodes exactly 2^N symbols, where N is the number of bits
per symbol in A.

This implies that every bitpattern of N bits is a valid encoding,
and every N-bit encoding decodes to a unique `BioSymbol`.

It defaults to `false`. Defining `iscomplete(::A) = true` for a new
`A <: Alphabet` enables certain optimisations.
"""
iscomplete(::Alphabet) = false

## Encoders & Decoders
"""
    abstract type AlphabetCode end

Trait for text parsing of biosymbols.
Currently, only `AsciiAlphabet` and `GenericAlphabetCode` is defined.
Define `AlphabetCode(::Alphabet)` to opt into non-generic fast paths
for parsing biosequences.

By default, `AlphabetCode(::AlphabetCode) === GenericAlphabetCode()`.

See also: [`ASCIIAlphabet`](@ref)
"""
abstract type AlphabetCode end

"""
    ASCIIAlphabet <: AlphabetCode

Singleton trait struct used to signify that all symbols of an `Alphabet`
are encoded from, and decoded to, an ASCII character, when parsing from,
or printing to a string.
For example, `DNAAlphabet{2}` is an `ASCIIAlphabet`, since all symbols in
this alphabet have a textual representation as ASCII.
"""
struct ASCIIAlphabet <: AlphabetCode end

"""
    GenericAlphabetCode <: AlphabetCode

The generic, fallback trait object used when parsing a `BioSequence` from,
or converting to, a string format.

By default, `AlphabetCode(::AlphabetCode) === GenericAlphabetCode()`.
"""
struct GenericAlphabetCode <: AlphabetCode end

AlphabetCode(::Alphabet) = GenericAlphabetCode()

"""
    abstract type EncodingScheme

Trait object used to guide dispatch when encoding data into `BioSequence`s.
When encoding data from a sequence of type `S`, into a sequence with alphabet
`A`, `EncodingScheme(::A, S)` should return some singleton struct `E <: EncodingScheme`,
which then guides dispatch to the most efficient method.
"""
abstract type EncodingScheme end

"""
    GenericEncoding <: EncodingScheme

Generic, fallback encoding scheme for arbitrary alphabets and iterables.

Generic methods should specialize on `::EncodingScheme`, not on `::GenericEncoding`,
to handle future new encoding schemes.
"""
struct GenericEncoding <: EncodingScheme end

EncodingScheme(::Alphabet, ::Type) = GenericEncoding()

"""
    CopyableEncoding <: EncodingScheme

This encoding scheme is used when copying data between `BioSequence`s of the same
`Alphabet`, or between alphabets where the binary encoding is compatible.
Note that this trait only implies compatible encoding for each individual biosymbol;
`BioSequence` subtypes may differ in how encodings are stored.
"""
struct CopyableEncoding <: EncodingScheme end

EncodingScheme(::A, ::Type{<:BioSequence{A}}) where {A <: Alphabet} = CopyableEncoding()

"""
    ASCIIEncoding <: EncodingScheme

Encoding scheme used when constructing a `BioSequence` with a `ASCIIAlphabet`
from an ASCII-compatible source.
For example, when constructing a DNA sequence from a `String`, the string
may be interpreted as ASCII codeunits, circumventing unicode decoding
altogether.
"""
struct ASCIIEncoding <: EncodingScheme end

function EncodingScheme(A::Alphabet, T::Type{<:AbstractString})
    return EncodingScheme(A, AlphabetCode(A), T)
end

function EncodingScheme(::Alphabet, ::GenericAlphabetCode, ::Type)
    return GenericEncoding()
end

function EncodingScheme(::Alphabet, ::ASCIIAlphabet, ::Type{<:Union{String, SubString{String}}})
    return ASCIIEncoding()
end

@static if isdefined(Base, :StringView)
    function EncodingScheme(
            ::Alphabet,
            ::ASCIIAlphabet,
            ::Type{<:Union{StringView, SubString{<:StringView}}}
        )
        ASCIIEncoding()
    end
end

"""
    EncodeError(A::Alphabet, val::T, idx::Integer)

Exception thrown when a `BioSymbol` cannot be encoded to a given [`Alphabet`](@ref).

# Examples
```
julia> try
           BioSequences.encode(DNAAlphabet{2}(), DNA_N)
       catch err
           println(err isa BioSequences.EncodeError)
       end
true
```
"""
struct EncodeError{A <: Alphabet, T} <: Exception
    val::T

    # One-based index of the source sequence where the bad symbol was found.
    # If not applicable (e.g. a single symbol was encoded), set to zero
    idx::Int
end

function EncodeError(::A, val::T, idx::Integer) where {A, T}
    return EncodeError{A, T}(val, Int(idx)::Int)
end

# TODO: Improve this
function Base.showerror(io::IO, err::EncodeError{A}) where {A}
    val = err.val
    char_repr = if val isa Integer && val < 0x80
        repr(val) * " (Char '" * Char(val) * "')"
    elseif val isa Union{AbstractString, AbstractChar}
        repr(val)
    else
        string(err.val)
    end
    return print(io, "cannot encode " * char_repr * " in ", A)
end

"""
    encode(A::Alphabet, s::BioSymbol)::UInt64

Encode `s` into the alphabet, or throw an `EncodeError(A, s)` if the `BioSymbol`
is not valid in the alphabet.
The encoding must be a `UInt64`, where the all but the lower N bits are unset,
where N is the number of bits per symbol.

Users should generally override [`tryencode`](@ref) instead.
"""
function encode(A::Alphabet, s::BioSymbol)
    y = tryencode(A, s)
    return y === nothing ? throw(EncodeError(A, s, 0)) : y
end

"""
    tryencode(::Alphabet, s::Biosymbol)::Union{Nothing, UInt64}

Encode `s` into the alphabet, or return `nothing` if the `BioSymbol`
is not valid in the alphabet.
The encoding must be a `UInt64`, where the all but the lower N bits are unset,
where N is the number of bits per symbol.
"""
function tryencode end

"""
	try_ascii_encode(::Alphabet, b::UInt8)::Union{Nothing, UInt64}

Encode the byte `b`, representing an ASCII character to the internal alphabet encoding.
This method is only used if `AlphabetCode(::A) == ASCIIAlphabet()`.

This function must return `nothing` if `b > 0x7f`, or if the ASCII character
with the codepoint `b` cannot be encoded into the alphabet.

See also: [`BioSequences.ascii_decode`](@ref), [`BioSequences.ASCIIAlphabet`](@ref)

```jldoctest
julia> try_ascii_encode(DNAAlphabet{2}(), 0x62)
0x0000000000000001

julia> try_ascii_encode(DNAAlphabet{2}(), UInt8('W'))

julia> try_ascii_encode(AminoAcidAlphabet(), UInt8('W'))
0x0000000000000011

julia> try_ascii_encode(DNAAlphabet{2}(), 0xf0)
```
"""
function try_ascii_encode end

"""
	ascii_decode(::Alphabet, enc::UInt64)::UInt8

Decode the internal encoding `enc` to the ASCII codepoint for that `BioSymbol`.
For example, `DNA_G` is encoded in `DNAAlphabet{4}()` as `UInt64(4)`, and so
`ascii_decode(DNAAlphabet{4}(), UInt64(4)) == UInt8('G')`.
This method is only used if `AlphabetCode(::A) == ASCIIAlphabet()`

The encoding can be assumed to be valid. Passing an invalid encoding
to this function may result in undefined behaviour.

See also: [`BioSequences.try_ascii_encode`](@ref), [`BioSequences.ASCIIAlphabet`](@ref)
"""
function ascii_decode end

"""
    decode(::Alphabet, x::UInt64)::BioSymbol

Decode internal representation (encoding) of the alphabet to the corresponding `BioSymbol`.

The encoding can be assumed to be valid. Passing an invalid encoding
to this function may result in undefined behaviour.

See also: [`BioSequences.encode`](@ref), [`Alphabet`](@ref)
"""
function decode end

"""
    to_decoding_source(scheme::EncodingScheme, x::T)::S

Construct an object suitable to acting as a decoding source in `scheme`.

This function is called on `x`, when constructing `BioSequence`s from `x`.
By default, it returns the input argument.
However, strings, and some certain `BioSequence` types can be cast as a
different, more easily parseable source without runtime cost.
"""
to_decoding_source(::EncodingScheme, x) = x

to_decoding_source(::ASCIIEncoding, x::Union{String, SubString{String}}) = codeunits(x)

@static if isdefined(Base, :StringView)
    function to_decoding_source(
            ::ASCIIEncoding,
            x::Union{StringView, SubString{<:StringView}},
        )
        codeunits(x)
    end
end

# TODO: InlineStrings and StringViews extension here.

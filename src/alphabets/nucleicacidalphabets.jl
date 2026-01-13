"""
    abstract type NucleicAcidAlphabet <: Alphabet end

Supertype of nucleic acids alphabets.

See also: [`DNAAlphabet`](@ref), [`RNAAlphabet`](@ref)
"""
abstract type NucleicAcidAlphabet <: Alphabet end

"""
    DNAAlphabet{N} <: NucleicAcidAlphabet

Dideoxyribonucleic acid alphabet. The parameter `N` must be an `Int`
and may be 2 or 4.
If N is 2, the alphabet can encode only A, C, G and T.
If N is 4, the alphabet can encode all instances of `DNA`.

The encoding of this alphabet is stable API. The encoding matches that
of `RNAAlphabet`.

### DNAAlphabet{2} encoding:
The encoding is a `UInt64` with the two lower bits given by the table below,
and all other bits zeroed

| Symbol  | Encoding |
| `DNA_A` | `00`     |
| `DNA_C` | `01`     |
| `DNA_G` | `10`     |
| `DNA_T` | `11`     |

### DNAAlphabet{4} encoding:
The encoding is a `UInt64` with the four lower bits given by the table
below, and all other bits zeroed.
This is guaranteed to be the same as `reinterpret(UInt8, symbol) % UInt64`.

| Symbol    | Encoding |
| `DNA_Gap` | `0000`   |
| `DNA_A`   | `0001`   |
| `DNA_C`   | `0010`   |
| `DNA_M`   | `0011`   |
| `DNA_G`   | `0100`   |
| `DNA_R`   | `0101`   |
| `DNA_S`   | `0110`   |
| `DNA_V`   | `0111`   |
| `DNA_T`   | `1000`   |
| `DNA_W`   | `1001`   |
| `DNA_Y`   | `1010`   |
| `DNA_H`   | `1011`   |
| `DNA_K`   | `1100`   |
| `DNA_D`   | `1101`   |
| `DNA_B`   | `1110`   |
| `DNA_N`   | `1111`   |

See also: [`RNAAlphabet`](@ref)
"""
struct DNAAlphabet{N} <: NucleicAcidAlphabet
    function DNAAlphabet{N}() where {N}
        if !(N isa Int) || (N != 2 && N != 4)
            error("In DNAAlphabet, N must be an Int equal to 2 or 4")
        end
        return new{N}()
    end
end

"""
    RNAAlphabet{N} <: NucleicAcidAlphabet

Ribonucleic acid alphabet. The parameter `N` must be an `Int`
and may be 2 or 4.
If N is 2, the alphabet can encode only A, C, G and U.
If N is 4, the alphabet can encode all instances of `RNA`.

The encoding of this alphabet is stable API. The encoding
matches that of `DNAAlphabet`.

### RNAAlphabet{2} encoding:
The encoding is a `UInt64` with the two lower bits given by the table below,
and all other bits zeroed

| Symbol  | Encoding |
| `RNA_A` | `00`     |
| `RNA_C` | `01`     |
| `RNA_G` | `10`     |
| `RNA_U` | `11`     |

### RNAAlphabet{4} encoding:
The encoding is a `UInt64` with the four lower bits given by the table
below, and all other bits zeroed.
This is guaranteed to be the same as `reinterpret(UInt8, symbol) % UInt64`.

| Symbol    | Encoding |
| `RNA_Gap` | `0000`   |
| `RNA_A`   | `0001`   |
| `RNA_C`   | `0010`   |
| `RNA_M`   | `0011`   |
| `RNA_G`   | `0100`   |
| `RNA_R`   | `0101`   |
| `RNA_S`   | `0110`   |
| `RNA_V`   | `0111`   |
| `RNA_U`   | `1000`   |
| `RNA_W`   | `1001`   |
| `RNA_Y`   | `1010`   |
| `RNA_H`   | `1011`   |
| `RNA_K`   | `1100`   |
| `RNA_D`   | `1101`   |
| `RNA_B`   | `1110`   |
| `RNA_N`   | `1111`   |

See also: [`DNAAlphabet`](@ref)
"""
struct RNAAlphabet{N} <: NucleicAcidAlphabet
    function RNAAlphabet{N}() where {N}
        if !(N isa Int) || (N != 2 && N != 4)
            error("In RNAAlphabet, N must be an Int equal to 2 or 4")
        end
        return new{N}()
    end
end

const TwoBitNuc = Union{DNAAlphabet{2}, RNAAlphabet{2}}
const FourBitNuc = Union{DNAAlphabet{4}, RNAAlphabet{4}}

Base.eltype(::Type{<:DNAAlphabet}) = DNA
Base.eltype(::Type{<:RNAAlphabet}) = RNA

Base.length(::DNAAlphabet{2}) = 4
Base.length(::DNAAlphabet{4}) = 16
Base.length(::RNAAlphabet{2}) = 4
Base.length(::RNAAlphabet{4}) = 16

BitsPerSymbol(::TwoBitNuc) = BitsPerSymbol{2}()
BitsPerSymbol(::FourBitNuc) = BitsPerSymbol{4}()

function Base.iterate(::A, state::UInt8 = 0x01) where {A <: TwoBitNuc}
    state > 0x08 && return nothing
    return (reinterpret(eltype(A), state), state << 1)
end

function Base.iterate(::A, state::UInt8 = 0x00) where {A <: FourBitNuc}
    state > 0x0f && return nothing
    return (reinterpret(eltype(A), state), state + 0x01)
end

function tryencode(::DNAAlphabet{2}, x::DNA)
    u = reinterpret(UInt8, x)
    count_ones(u) == 1 || return nothing
    return trailing_zeros(u) % UInt64
end

function tryencode(::RNAAlphabet{2}, x::RNA)
    u = reinterpret(UInt8, x)
    count_ones(u) == 1 || return nothing
    return trailing_zeros(u) % UInt64
end

tryencode(::DNAAlphabet{4}, x::DNA) = reinterpret(UInt8, x) % UInt64
tryencode(::RNAAlphabet{4}, x::RNA) = reinterpret(UInt8, x) % UInt64

function decode(::A, x::UInt64) where {A <: TwoBitNuc}
    return reinterpret(eltype(A), 0x01 << (x & UInt64(3)))
end

function decode(::A, x::UInt64) where {A <: FourBitNuc}
    return reinterpret(eltype(A), x % UInt8)
end

iscomplete(::Union{TwoBitNuc, FourBitNuc}) = true

AlphabetCode(::Union{DNAAlphabet, RNAAlphabet}) = ASCIIAlphabet()

struct TwoToFour <: EncodingScheme end

struct FourToTwo <: EncodingScheme end

const (
    DNA2_ASCII_ENCODE_LUT,
    RNA2_ASCII_ENCODE_LUT,
    DNA4_ASCII_ENCODE_LUT,
    RNA4_ASCII_ENCODE_LUT,
) = let
    vectors = [fill(0xff, 256) for i in 1:4]
    for (vector, A) in zip(vectors, Any[DNAAlphabet{2}(), RNAAlphabet{2}(), DNAAlphabet{4}(), RNAAlphabet{4}()])
        for symbol in A
            encoding = UInt8(encode(A, symbol))
            for char in [lowercase(Char(symbol)), uppercase(Char(symbol))]
                vector[Int(char) + 1] = encoding
            end
        end
    end
    map(Tuple, vectors)
end

function try_ascii_encode(::DNAAlphabet{2}, x::UInt8)
    encoding = @inbounds DNA2_ASCII_ENCODE_LUT[x + 1]
    return encoding == 0xff ? nothing : (encoding % UInt64)
end

function try_ascii_encode(::RNAAlphabet{2}, x::UInt8)
    encoding = @inbounds RNA2_ASCII_ENCODE_LUT[x + 1]
    return encoding == 0xff ? nothing : (encoding % UInt64)
end

function try_ascii_encode(::DNAAlphabet{4}, x::UInt8)
    encoding = @inbounds DNA4_ASCII_ENCODE_LUT[x + 1]
    return encoding == 0xff ? nothing : (encoding % UInt64)
end

function try_ascii_encode(::RNAAlphabet{4}, x::UInt8)
    encoding = @inbounds RNA4_ASCII_ENCODE_LUT[x + 1]
    return encoding == 0xff ? nothing : (encoding % UInt64)
end

const (
    DNA2_ASCII_DECODE_LUT,
    RNA2_ASCII_DECODE_LUT,
    DNA4_ASCII_DECODE_LUT,
    RNA4_ASCII_DECODE_LUT,
) = let
    alphabets = Any[DNAAlphabet{2}(), RNAAlphabet{2}(), DNAAlphabet{4}(), RNAAlphabet{4}()]
    vectors = map(A -> fill(0xff, length(A)), alphabets)
    for (A, vector) in zip(alphabets, vectors)
        for symbol in A
            encoding = UInt8(encode(A, symbol))
            vector[encoding + 1] = Char(symbol)
        end
    end
    map(Tuple, vectors)
end

function ascii_decode(::DNAAlphabet{2}, x::UInt64)
    return @inbounds DNA2_ASCII_DECODE_LUT[(x % Int) + 1]
end

function ascii_decode(::RNAAlphabet{2}, x::UInt64)
    return @inbounds RNA2_ASCII_DECODE_LUT[(x % Int) + 1]
end

function ascii_decode(::DNAAlphabet{4}, x::UInt64)
    return @inbounds DNA4_ASCII_DECODE_LUT[(x % Int) + 1]
end

function ascii_decode(::RNAAlphabet{4}, x::UInt64)
    return @inbounds RNA4_ASCII_DECODE_LUT[(x % Int) + 1]
end

"""
    AminoAcidAlphabet <: Alphabet

Alphabet that encompasses all instances of `AminoAcid`, and coded
in 8 bits.

See also: [`Alphabet`](@ref)
"""
struct AminoAcidAlphabet <: Alphabet end

BitsPerSymbol(::AminoAcidAlphabet) = BitsPerSymbol{8}()

Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid

Base.length(::AminoAcidAlphabet) = 28

function Base.iterate(::AminoAcidAlphabet, state::UInt8=0x00)
    state > 0x1b && return nothing
    (reinterpret(AminoAcid, state), state + 0x01)
end

function tryencode(::AminoAcidAlphabet, x::AminoAcid)
    reinterpret(UInt8, x) % UInt64
end

decode(::AminoAcidAlphabet, x::UInt64) = reinterpret(AminoAcid, x % UInt8)

AlphabetCode(::AminoAcidAlphabet) = ASCIIAlphabet()

const (AA_ASCII_ENCODE_LUT, AA_ASCII_DECODE_LUT) = let
    encode_lut = fill(0xff, 256)
    decode_lut = fill(0xff, length(AminoAcidAlphabet()))

    for symbol in AminoAcidAlphabet()
        decode_lut[reinterpret(UInt8, symbol) + 1] = UInt8(Char(symbol))
        for char in [lowercase(Char(symbol)), uppercase(Char(symbol))]
            encode_lut[UInt8(char) + 1] = reinterpret(UInt8, symbol)
        end
    end

    (Tuple(encode_lut), Tuple(decode_lut))
end

function try_ascii_encode(::AminoAcidAlphabet, x::UInt8)
    encoding = @inbounds AA_ASCII_ENCODE_LUT[x + 1]
    encoding == 0xff ? nothing : (encoding % UInt64) 
end

function ascii_decode(::AminoAcidAlphabet, x::UInt64)
    @inbounds AA_ASCII_DECODE_LUT[(x % Int) + 1]
end
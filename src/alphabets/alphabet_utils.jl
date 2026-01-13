# This contains 128 UInt8, representing a bitmask of which alphabets can encode each
# ASCII codepoint.
# So, GUESS_ALPHABET_LUT[Int('C') + 1] is a bitmask containing the alphabets
# that can contain 'C' (currently: all alphabets)
# The bits from lowest to highest are given in the array below
const GUESS_ALPHABET_LUT = let
    v = zeros(UInt8, 128)
    for (offset, A) in [
            (0, DNAAlphabet{2}()),
            (1, RNAAlphabet{2}()),
            (2, DNAAlphabet{4}()),
            (3, RNAAlphabet{4}()),
            (4, AminoAcidAlphabet()),
        ]
        for symbol in A
            for byte in [UInt8(uppercase(Char(symbol))), UInt8(lowercase(Char(symbol)))]
                v[byte + 1] |= 0x01 << offset
            end
        end
    end
    Tuple(v)
end

"""
    guess_alphabet(s::Union{AbstractString, AbstractVector{UInt8}}) -> Union{Integer, Alphabet}

Pick an `Alphabet` that can encode input `s`.  If no `Alphabet` can, return the index of the first
byte of the input which is not encodable in any alphabet.
This function only knows about the alphabets listed below. If multiple alphabets are possible,
pick the first from the order below (i.e. `DNAAlphabet{2}()` if possible, otherwise `RNAAlphabet{2}()` etc).
1. `DNAAlphabet{2}()`
2. `RNAAlphabet{2}()`
3. `DNAAlphabet{4}()`
4. `RNAAlphabet{4}()`
5. `AminoAcidAlphabet()`

!!! warning
    The functions `bioseq` and `guess_alphabet` are intended for use in interactive
    sessions, and are not suitable for use in packages or non-ephemeral work.
    They are type unstable, and their heuristics **are subject to change** in minor versions.

# Examples
```jldoctest
julia> guess_alphabet("AGGCA")
DNAAlphabet{2}()

julia> guess_alphabet("WKLQSTV")
AminoAcidAlphabet()

julia> guess_alphabet("QAWT+!")
5

julia> guess_alphabet("UAGCSKMU")
RNAAlphabet{4}()
```
"""
function guess_alphabet(v::AbstractVector{UInt8})
    possibilities = 0x1f
    for (index, byte) in pairs(v)
        lut_byte = @inbounds GUESS_ALPHABET_LUT[(byte & 0x7f) + 0x01]
        possibilities &= (lut_byte * (byte < 0x80))
        iszero(possibilities) && return index
    end
    @assert !iszero(possibilities) # We checked that in the loop above
    if !iszero(possibilities &     0b00001)
        DNAAlphabet{2}()
    elseif !iszero(possibilities & 0b00010)
        RNAAlphabet{2}()
    elseif !iszero(possibilities & 0b00100)
        RNAAlphabet{2}()
    elseif !iszero(possibilities & 0b01000)
        RNAAlphabet{2}()
    else
        AminoAcidAlphabet()
    end
end

guess_alphabet(s::AbstractString) = guess_alphabet(codeunits(s))

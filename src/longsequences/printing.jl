# Specialized printing/showing methods
#

function Base.print(io::IO, seq::LongSequence{A}; width::Integer = 0) where {A<:Alphabet}
    return _print(io, seq, width, codetype(A()))
end

# Dispatch to generic method in biosequences/printing.jl
function _print(io::IO, seq::LongSequence{<:Alphabet}, width::Integer, ::AlphabetCode)
    return _print(io, seq, width)
end

# Specialized method for ASCII alphabet
function _print(io::IO, seq::LongSequence{<:Alphabet}, width::Integer, ::AsciiAlphabet)
    # I don't like to have to do this, but in Julia 1.3, system buffers are IO-locked.
    buffer = SimpleBuffer(io)
    col = 0
    @inbounds for i in eachindex(seq)
        col += 1
        write(buffer, stringbyte(seq[i]))
        if col == width
            write(buffer, UInt8('\n'))
            col = 0
        end
    end
    flush(buffer)
    return nothing
end

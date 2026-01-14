function write_into!(v::AbstractVector{UInt8}, seq::BioSequence)
    Base.require_one_based_indexing(v)
    @boundscheck if length(v) != length(seq)
        throw(DimensionMismatch("Length of vector must match length of biosequence"))
    end
    for i in eachindex(seq)
        enc = unsafe_get_encoding(seq, i)
        @inbounds v[i] = ascii_decode(Alphabet(seq), enc)
    end
    return v
end

Base.summary(seq::BioSequence{<:DNAAlphabet}) = string(length(seq), "nt ", "DNA Sequence")
Base.summary(seq::BioSequence{<:RNAAlphabet}) = string(length(seq), "nt ", "RNA Sequence")
Base.summary(seq::BioSequence{<:AminoAcidAlphabet}) = string(length(seq), "aa ", "Amino Acid Sequence")

function Base.String(seq::BioSequence)
    return biosequence_to_string(String, AlphabetCode(Alphabet(seq)), seq)
end

function biosequence_to_string(::Type{String}, ::ASCIIAlphabet, seq::BioSequence)
    v = write_into!(Vector{UInt8}(undef, length(seq)), seq)
    return String(v)
end

function biosequence_to_string(::Type{String}, ::AlphabetCode, seq::BioSequence)
    io = IOBuffer()
    for symbol in seq
        write(io, symbol)
    end
    return takestring!(io)
end

function Base.print(io::IO, seq::BioSequence)
    return _print(io, AlphabetCode(Alphabet(seq)), seq)
end

function _print(io::IO, ::AlphabetCode, seq::BioSequence)
    for symbol in seq
        print(io, symbol)
    end
    return
end

function _print(io::IO, ::ASCIIAlphabet, seq::BioSequence)
    chunk_size = 1024
    buffer = Memory{UInt8}(undef, min(length(seq), chunk_size))
    for i in 1:chunk_size:max(1, (length(seq) - chunk_size + 1))
        seq_view = @inbounds view(seq, i:min(i + chunk_size - 1, length(seq)))
        buffer_view = @inbounds view(buffer, 1:length(seq_view))
        write_into!(buffer_view, seq_view)
        write(io, buffer_view)
    end
    return
end

Base.show(io::IO, seq::BioSequence) = showcompact(io, seq)

function Base.show(io::IO, ::MIME"text/plain", seq::BioSequence)
    println(io, summary(seq), ':')
    return showcompact(io, seq)
end

function showcompact(io::IO, seq::BioSequence)
    # don't show more than this many characters
    # to avoid filling the screen with junk
    return if isempty(seq)
        print(io, "< EMPTY SEQUENCE >")
    else
        width = displaysize()[2]
        if length(seq) > width
            half = div(width, 2)
            for i in 1:(half - 1)
                print(io, seq[i])
            end
            print(io, '…')
            for i in (lastindex(seq) - half + 2):lastindex(seq)
                print(io, seq[i])
            end
        else
            print(io, seq)
        end
    end
end

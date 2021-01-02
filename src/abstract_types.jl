#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

"""
    abstract type FiniteField end
"""
abstract type FiniteField end

"""
    abstract type AbstractCode end
"""
abstract type AbstractCode end

"""
    NonStaticAbstractWord{N, T} = Union{NTuple{N, T}, AbstractVector{T}} where {N, T}
"""
NonStaticAbstractWord{N, T} = Union{NTuple{N, T}, AbstractVector{T}} where {N, T}

"""
    mutable struct Word{N, T}
    Word(w::NTuple{N, T})
    Word(w::AbstractVector{T})
    Word(i::T...)

A Word is an `StaticArrays.MVector` which is efficient (like tuple) yet mutable.
"""
mutable struct Word{N, T}
    w::MVector{N, T}
    
    Word(w::NTuple{N, T}) where {N, T} = new{N, T}(MVector{N, T}(w))
    Word(w::AbstractVector{T}) where {T} =
        (len = length(w); new{len, T}(MVector{len, T}(w)))
    Word(i::T...) where {T} =
        (len = length(i); new{len, T}(MVector{len, T}(i)))
end

# Indexing Interface
Base.getindex(W::Word{N, T}, i::Int) where {N, T} = getindex(W.w, i)
Base.setindex!(W::Word{N, T}, v, i::Int) where {N, T} = setindex!(W.w, v, i)
Base.firstindex(W::Word{N, T}) where {N, T} = firstindex(W.w)
Base.lastindex(W::Word{N, T}) where {N, T} = lastindex(W.w)

# Abstract Array Interface
Base.size(::Word{N, T}) where {N, T} = tuple(N)
Base.length(::Word{N, T}) where {N, T} = N

"""
    AbstractWord{N} = Union{NonStaticAbstractWord{N, T}, MVector{T}} where {T}
"""
AbstractWord{N, T} = Union{Word{N, T}, NonStaticAbstractWord{N, T}, MVector{N, T}} where {N, T}

isword(w) = w isa Word
isword(i...) = isword(i)
isabstractword(w) = w isa AbstractWord
isabstractword(i...) = isabstractword(i)

"""
    struct Alphabet <: AbstractCode
        
Has the parameter Σ, which is the alphabet; a collection of strings, characters, symbols, or Ints.

---

    Alphabet(Σ::AbstractArray)

A constructor method for the struct Alphabet.  Takes an array of letters in the alphabet, and attempts to parse them as 64-bit Ints.
    
---

    Alphabet(Σ::AbstractArray)
	Alphabet(Σ::AbstractString)

A constructor method for the struct Alphabet.  Takes in a symbols and splits it into constituent characters.  Those symbols are the letters in the alphabet.  Will attempt to parse these as 64-bit Ints.
"""
struct Alphabet <: AbstractVector{Symbol}
    Σ::AbstractVector{Symbol}

    Alphabet(Σ::Union{Vector{T}, String}) where {T} = new(ensure_symbolic(unique(Σ)))
end # end struct

# Indexing Interface
Base.getindex(A::Alphabet, i::Int) = getindex(A.Σ, i)
Base.setindex!(A::Alphabet, v, i::Int) = setindex!(A.Σ, v, i)
Base.firstindex(A::Alphabet) = firstindex(A.Σ)
Base.lastindex(A::Alphabet) = lastindex(A.Σ)

# Abstract Array Interface
Base.size(A::Alphabet) = size(A.Σ)
Base.length(A::Alphabet) = length(A.Σ)
Base.rand(A::Alphabet) = rand(A.Σ)

"""
    gensym(q::Int) -> Vector{Symbol}

Generates a vector of unique symbols of length q.
"""
Base.gensym(q::Int) = Symbol[gensym() for _ in 1:q]

"""
    genalphabet(q::Int)

Generates an alphabet of q unique symbols.
"""
genalphabet(q::Int) = Alphabet(gensym(q))

"""
    struct UniverseParameters <: AbstractCode
    
Defines a structure for the messages in the code.  Parameters are the alphabet `Σ`, size of alphabet `q`, and block length `n`

    UniverseParameters(Σ::Alphabet, n::Int)
    UniverseParameters(Σ::AbstractArray, n::Int)
    UniverseParameters(Σ::Alphabet, q::Int, n::Int)
    UniverseParameters(Σ::AbstractArray, q::Int, n::Int)
    UniverseParameters(q::Int, n::Int)
    
An inner constructor function on the structure `UniverseParameters`.
"""
struct UniverseParameters <: AbstractCode
    Σ::Alphabet
    q::Int # size of alphabet
    n::Int # block length
    
    UniverseParameters(Σ::Alphabet, n::Int) = new(Σ, length(q), n)
    UniverseParameters(Σ::Alphabet, q::Int, n::Int) = new(Σ, q, n)
    UniverseParameters(Σ::AbstractVector{T}, n::Int) where {T} = new(Alphabet(Σ), length(Σ), n)
    UniverseParameters(Σ::AbstractVector{T}, q::Int, n::Int) where {T} = new(Σ, q, n)
    UniverseParameters(q::Int, n::Int) = new(genalphabet(q), q, n)
end

# Iteration inteface functions
Base.length(iter::UniverseParameters) = big(iter.q)^iter.n
Base.eltype(iter::UniverseParameters) = NTuple{iter.n, Symbol}

# Other Base interface functions for UniverseParameters
Base.rand(𝒰::UniverseParameters) = ntuple(_ -> rand(𝒰.Σ), 𝒰.n)

"""
	rand(𝒰::UniverseParameters, C::AbstractArray) -> Tuple

Given universe parameters 𝒰 and a code C, return a tuple including
  - A random word in C;
  - A random letter in the alphabet; and
  - A random index in the block length.
"""
Base.rand(𝒰::UniverseParameters, C::AbstractArray) = tuple(Word(rand(C)), rand(𝒰.Σ), rand(1:𝒰.n))

"""
    struct CodeUniverseIterator <: AbstractCode

A structure used to iterate through a code universe, with specified universe parameters.  E.g.,

    for c in CodeUniverseIterator(["a", "b", "c"], 3)
        println(c)
    end

Fields:
  - 𝒰::UniverseParameters

Methods:

    CodeUniverseIterator(𝒰::UniverseParameters)
    CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, q::Int, n::Int)
    CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, n::Int)
    CodeUniverseIterator(q::Int, n::Int)

"""
struct CodeUniverseIterator <: AbstractCode
    𝒰::UniverseParameters

    function CodeUniverseIterator(𝒰::UniverseParameters)
        Σ = 𝒰.Σ
        return Iterators.product(Vector{eltype(Σ)}[Σ for _ in 1:𝒰.n]...)
    end

    CodeUniverseIterator(Σ::Union{Alphabet, AbstractVector{T}}, q::Int, n::Int) where {T} =
        CodeUniverseIterator(UniverseParameters(Σ, q, n))
    CodeUniverseIterator(Σ::Union{Alphabet, AbstractVector{T}}, n::Int) where {T} =
        CodeUniverseIterator(UniverseParameters(Σ, n))
    CodeUniverseIterator(q::Int, n::Int) =
        CodeUniverseIterator(UniverseParameters(q, n))
end

"""
    struct CodeUniverse <: AbstractCode
    
Defines a structure for the messages in the code.  Parameters are the abstract array of messages `𝒰`, and the length of the messages `n`.

    CodeUniverse(𝒰::AbstractArray, Σ::Alphabet)
    
An inner constructor function on the structure `CodeUniverse`.
"""
struct CodeUniverse <: AbstractCode
    𝒰::AbstractVector{Word{N, T}} where {N, T}
    Σ::Alphabet
    q::Int # alphabet size
    n::Int # block length
    
    function CodeUniverse(𝒰::AbstractVector{Word{N, T}}, Σ::Alphabet) where {N, T}
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(𝒰) || throw(error(message_length_error))
    
        new(𝒰, Σ, length(Σ), length(first(𝒰)))
    end # end constructor function
end

"""
    struct Rounding end

An abstract structure used as a parameter in a non-rounding method of `hamming_bound`.  Use `no_round` for this; a predefiend variable of the structure `Rounding`.
"""
struct Rounding end
no_round = Rounding()

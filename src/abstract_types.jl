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
    struct Alphabet <: AbstractCode
        
Has the parameter Σ, which is the alphabet; a collection of strings, characters, symbols, or integers.

---

    Alphabet(Σ::AbstractArray)

A constructor method for the struct Alphabet.  Takes an array of letters in the alphabet, and attempts to parse them as 64-bit integers.
    
---

    Alphabet(Σ::AbstractArray)
	Alphabet(Σ::AbstractString)

A constructor method for the struct Alphabet.  Takes in a symbols and splits it into constituent characters.  Those symbols are the letters in the alphabet.  Will attempt to parse these as 64-bit integers.
"""
struct Alphabet <: AbstractVector{Symbol}
    Σ::AbstractVector{Symbol}
    
    function Alphabet(Σ::AbstractArray)
        Σ = ensure_symbolic(unique(Σ))
        
        new(Σ)
    end # end constructor function

    function Alphabet(Σ::AbstractString)
        Σ = ensure_symbolic(collect(unique(Σ)))
        
        new(Σ)
    end # end constructor function
end # end struct

# Indexing Interface
Base.getindex(A::Alphabet, i::Integer) = getindex(A.Σ, i)
Base.setindex!(A::Alphabet, v, i::Integer) = setindex!(A.Σ, v, i)
Base.firstindex(A::Alphabet) = firstindex(A.Σ)
Base.lastindex(A::Alphabet) = lastindex(A.Σ)

# Abstract Array Interface
Base.size(A::Alphabet) = size(A.Σ)
Base.getindex(A::Alphabet, i::Integer) = getindex(A.Σ, i)
Base.setindex!(A::Alphabet, v, i::Integer) = setindex(A.Σ, v, i)
Base.rand(A::Alphabet) = rand(A.Σ)

"""
    struct UniverseParameters <: AbstractCode
    
Defines a structure for the messages in the code.  Parameters are the alphabet `Σ`, size of alphabet `q`, and block length `n`

    UniverseParameters(Σ::Alphabet, n::Integer)
    UniverseParameters(Σ::AbstractArray, n::Integer)
    UniverseParameters(Σ::Alphabet, q::Integer, n::Integer)
    UniverseParameters(Σ::AbstractArray, q::Integer, n::Integer)
    UniverseParameters(q::Integer, n::Integer)
    
An inner constructor function on the structure `UniverseParameters`.
"""
struct UniverseParameters <: AbstractCode
    Σ::Alphabet
    q::Integer
    n::Integer # block length
    
    function UniverseParameters(Σ::Alphabet, n::Integer)
        q = length(Σ)
        
        new(Σ, q, n)
    end
    
    UniverseParameters(Σ::Alphabet, q::Integer, n::Integer) = new(Σ, q, n)
    
    function UniverseParameters(Σ::AbstractArray, n::Integer)
        Σ = Alphabet(Σ)
        q = length(Σ)
        
        new(Σ, q, n)
    end
    
    UniverseParameters(Σ::AbstractArray, q::Integer, n::Integer) = new(Σ, q, n)
    
    function UniverseParameters(q::Integer, n::Integer)
        Σ = Alphabet([gensym() for i in 1:q])
        
        new(Σ, q, n)
    end
end

# Iteration inteface functions
Base.length(iter::UniverseParameters) = iter.q^iter.n
Base.eltype(iter::UniverseParameters) = NTuple{iter.n, Symbol}

# Other Base interface functions for UniverseParameters
Base.rand(𝒰::UniverseParameters) = rand(𝒰.Σ)

"""
	rand(𝒰::UniverseParameters, C::AbstractArray) -> Tuple

Given universe parameters 𝒰 and a code C, return a tuple including
  - A random word in C;
  - A random letter in the alphabet; and
  - A random index in the block length.
"""
Base.rand(𝒰::UniverseParameters, C::AbstractArray) = rand.((C, 𝒰.Σ, 1:𝒰.n))

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
    CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, q::Integer, n::Integer)
    CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, n::Integer)
    CodeUniverseIterator(q::Integer, n::Integer)

"""
struct CodeUniverseIterator <: AbstractCode
    𝒰::UniverseParameters

    function CodeUniverseIterator(𝒰::UniverseParameters)
        return reshape(collect(Iterators.product([𝒰.Σ for _ in 1:𝒰.n]...)), :)
    end

    CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, q::Integer, n::Integer) =
        CodeUniverseIterator(UniverseParameters(Σ, q, n))
    CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, n::Integer) =
        CodeUniverseIterator(UniverseParameters(Σ, n))
    CodeUniverseIterator(q::Integer, n::Integer) =
        CodeUniverseIterator(UniverseParameters(q, n))
end

"""
    struct CodeUniverse <: AbstractCode
    
Defines a structure for the messages in the code.  Parameters are the abstract array of messages `𝒰`, and the length of the messages `n`.

    CodeUniverse(𝒰::AbstractArray, Σ::Alphabet)
    
An inner constructor function on the structure `CodeUniverse`.
"""
struct CodeUniverse <: AbstractCode
    𝒰::AbstractArray
    Σ::Alphabet
    q::Integer
    n::Integer # block length
    
    
    function CodeUniverse(𝒰::AbstractArray, Σ::Alphabet)
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(𝒰) || throw(error("$(message_length_error)"))
    
        q = length(Σ)
        n = length(rand(𝒰)) # choose arbitrary message in the list of messages
    
        new(𝒰, Σ, q, n)
    end # end constructor function
end

"""
    struct Rounding end

An abstract structure used as a parameter in a non-rounding method of `hamming_bound`.  Use `no_round` for this; a predefiend variable of the structure `Rounding`.
"""
struct Rounding end
no_round = Rounding()

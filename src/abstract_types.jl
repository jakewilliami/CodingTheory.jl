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

A constructor method for the struct Alphabet.  Takes in a symbols and splits it into constituent characters.  Those symbols are the letters in the alphabet.  Will attempt to parse these as 64-bit integers.
"""
struct Alphabet <: AbstractVector{Symbol}
    values::AbstractArray{Symbol}
    
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

Base.getindex(A::Alphabet, i) = A.values[i]
Base.setindex!(A::Alphabet, v, i) = (A.values[i] = v)
Base.firstindex(A::Alphabet) = firstindex(A.values)
Base.lastindex(A::Alphabet) = lastindex(A.values)

# Abstract Array Interface

Base.size(A::Alphabet) = size(A.values)
Base.getindex(A::Alphabet, i::Int) = getindex(A.values, i)
Base.setindex!(A::Alphabet, v, i::Int) = setindex(A.values, v, i)

"""
    struct CodeUniverse <: AbstractCode
    
Defines a structure for the messages in the code.  Parameters are the abstract array of messages `C`, and the length of the messages `n`.

    CodeUniverse(C::AbstractArray, Σ::Alphabet)
    
An inner constructor function on the structure `CodeUniverse`.
"""
struct CodeUniverse <: AbstractCode
    C::AbstractArray
    Σ::Alphabet
    q::Integer
    n::Integer # block length
    
    
    function CodeUniverse(C::AbstractArray, Σ::Alphabet)
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(C) || throw(error("$(message_length_error)"))
        
        q = length(Σ)
        n = length(rand(C)) # choose arbitrary message in the list of messages
        
        new(C, Σ, q, n)
    end # end constructor function
end

"""
    struct UniverseParameters <: AbstractCode
    
Defines a structure for the messages in the code.  Parameters are the alphabet `Σ`, size of alphabet `q`, and block length `n`

    UniverseParameters(Σ::Alphabet, n::Integer)
    UniverseParameters(Σ::AbstractArray, n::Integer)
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
    
    function UniverseParameters(Σ::AbstractArray, n::Integer)
        Σ = Alphabet(Σ)
        q = length(Σ)
        
        new(Σ, q, n)
    end
    
    function UniverseParameters(q::Integer, n::Integer)
        Σ = Alphabet([gensym() for i in 1:q])
        
        new(Σ, q, n)
    end
end

function Base.iterate(iter::UniverseParameters)
	if iszero(iter.n) || iszero(iter.q)
		return nothing
	end
	
	element = ntuple(_ -> first(iter.Σ), iter.n)
	
	return element, BigInt(0)
end

# "Induction" iterate method
function Base.iterate(iter::UniverseParameters, state::BigInt)
    state += 1
    if state >= iter.q ^ iter.n
        return nothing
    end
    word = ntuple(i -> iter.Σ[BigInt(floor(state / iter.q^(i-1))) % iter.q + 1], iter.n);
	return word, state
end

Base.length(iter::UniverseParameters) = iter.q^iter.n
Base.eltype(iter::UniverseParameters) = Tuple{Symbol}

"""
    struct Rounding end

An abstract structure used as a parameter in a non-rounding method of `hamming_bound`.  Use `no_round` for this; a predefiend variable of the structure `Rounding`.
"""
struct Rounding end
no_round = Rounding()

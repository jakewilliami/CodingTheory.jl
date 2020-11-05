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
struct Alphabet <: AbstractCode
    Σ::Union{AbstractString, AbstractChar, Symbol, Integer, AbstractArray{T}} where T
    
    function Alphabet(Σ::AbstractArray)
        Σ = ensure_symbolic(unique(Σ))
        
        new(Σ)
    end # end constructor function

    function Alphabet(Σ::AbstractString)
        Σ = ensure_symbolic(collect(unique(Σ)))
        
        new(Σ)
    end # end constructor function
end # end struct

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
        q = length(Σ.Σ)
        
        new(Σ, q, n)
    end
    
    function UniverseParameters(Σ::AbstractArray, n::Integer)
        Σ = Alphabet(Σ)
        q = length(Σ.Σ)
        
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
	
	element = ntuple(_ -> first(iter.Σ.Σ), iter.n)
	
	return element, element
end

function Base.iterate(iter::UniverseParameters, state)
	word = [i for i in state]
	i = 1
	
	while isequal(word[i], last(iter.Σ.Σ))
		word[i] = first(iter.Σ.Σ)
		
		if isequal(i, length(word))
			return nothing
		end
		
		i += 1
	end

	alphabet_index = findfirst(isequal(word[i]), iter.Σ.Σ)
	word[i] = iter.Σ.Σ[alphabet_index + 1]
	element = ntuple(x -> word[x], iter.n)
	
	return element, element
end

Base.length(iter::UniverseParameters) = iter.q^iter.n
Base.eltype(iter::UniverseParameters) = Tuple{Symbol}

"""
    struct Rounding end

An abstract structure used as a parameter in a non-rounding method of `hamming_bound`.  Use `no_round` for this; a predefiend variable of the structure `Rounding`.
"""
struct Rounding end
no_round = Rounding()

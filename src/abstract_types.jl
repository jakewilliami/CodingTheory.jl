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
        Σ = unique(Σ)
        
        try
            Σ = parse.(Int, Σ)
        catch
        end
        
        new(Σ)
    end # end constructor function

    function Alphabet(Σ::AbstractString)
        Σ = collect(unique(Σ))
        
        try
            Σ = parse.(Int, Σ)
        catch
            Σ = Symbol.(Σ)
        end
        
        new(Σ)
    end # end constructor function
end # end struct

"""
    struct Messages <: AbstractCode
    
Defines a structure for the messages in the code.  Parameters are the abstract array of messages `ℳ`, and the length of the messages `block_length`.

    Messages(ℳ::AbstractArray)
    
An inner constructor function on the structure `Messages`.  Ensures the messages are all of equal length.
"""
struct Messages <: AbstractCode
    ℳ::AbstractArray
    block_length::Integer
    
    function Messages(ℳ::AbstractArray)
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(ℳ) || throw(error("$(message_length_error)"))
        
        block_length = length(rand(ℳ)) # choose arbitrary message in the list of messages
        
        new(ℳ, block_length)
    end # end constructor function
end

"""
    struct Rounding end

An abstract structure used as a parameter in a non-rounding method of `hamming_bound`.  Use `no_round` for this; a predefiend variable of the structure `Rounding`.
"""
struct Rounding end
no_round = Rounding()

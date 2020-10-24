#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
abstract type FiniteField end
abstract type AbstractCode end
    
struct Alphabet <: AbstractCode
    Σ::Union{AbstractArray, AbstractString}
    
    function Alphabet(Σ::AbstractArray)
        try
            Σ = parse.(Int, Σ)
        catch
        end
        
        new(Σ)
    end # end constructor function
    
    function Alphabet(Σ::AbstractString)
        Σ = collect(Σ)
        
        try
            Σ = parse.(Int, Σ)
        catch
            Σ = string.(Σ)
        end
        
        new(Σ)
    end # end constructor function
end # end struct

struct Messages <: AbstractCode
    ℳ::AbstractArray
    length::Integer
    
    function Messages(ℳ::AbstractArray)
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(ℳ) || throw(error("$(message_length_error)"))
        
        length = length(ℳ[rand(1:length(ℳ))]) # choose arbitrary message in the list of messages
        
        new(ℳ, length)
    end # end constructor function
end

struct Rounding end
no_round = Rounding()

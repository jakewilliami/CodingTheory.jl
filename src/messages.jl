#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "utils.jl"))

struct Alphabet
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

struct Messages
    ℳ::AbstractArray
    length::Integer
    
    function Messages(ℳ::AbstractArray)
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(ℳ) || throw(error("$(message_length_error)"))
        
        length = length(ℳ[rand(1:length(ℳ))]) # choose arbitrary message in the list of messages
        
        new(ℳ, length)
    end # end constructor function
end

#=
n is the word length
q is the number of symbols in the code
M is the size/number of elements in the code
=#
function rate(q::Integer, M::Integer, n::Integer)
    return log(q, M) / n
end

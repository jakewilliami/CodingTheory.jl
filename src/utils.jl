#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

#=
Check that all elements in a list are of equal length.
=#
@inline function __allequal_length_(A::AbstractArray)::Bool
    length(A) < 2 && return true
	
    @inbounds for i in 2:length(A)
        isequal(length(A[1]), length(A[i])) || return false
    end
	
    return true
end

@inline function __allequal_length_(a...)::Bool
    A = [a...]
    return __allequal_length_(A)
end

#=
Check that all elements in a list are equal to each other.
=#
@inline function __allequal(A::AbstractArray)::Bool
    length(A) < 2 && return true
    
    @inbounds for i in 2:length(A)
        A[1] ≠ A[i] && return false
    end
    
    return true
end

@inline function __allequal(a...)::Bool
    A = [a...]
    return __allequal(A)
end

#=
Check that all elements in a list are distinct from every other element in the list.
=#
@inline function __aredistinct(A::AbstractArray)::Bool
    length(A) < 2 && return true
    
    while ! iszero(length(A))
        a = pop!(A)
        a ∉ A || return false
    end
    
    return true
end

@inline function __aredistinct(a...)::Bool
    A = [a...]
    return __aredistinct(A)
end

#=
Check that all elements in a list are less than a given x.
=#
@inline function __arelessthan(x::Number, A::AbstractArray)::Bool
    @inbounds for a in A
        a < x && return true
    end
    
    return false
end

@inline function __arelessthan(x::Number, a::Number...)::Bool
    A = [a...]
    return __arelessthan(x, A)
end

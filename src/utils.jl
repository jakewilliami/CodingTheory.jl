#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
	
using LinearAlgebra: I

function displaymatrix(M::AbstractArray)
    return show(IOContext(stdout, :limit => true, :compact => true, :short => true), "text/plain", M); print("\n")
end

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

#=
Takes in an array and a word.  As long as the word does not mean that the distance is smaller than d, we add w to the array
=#
function __push_if_allowed!(C::AbstractArray{T}, w::T, d::Integer) where T
	isempty(C) && push!(C, w)
	
	for c in C
		if hamming_distance(c, w) < d
			return nothing
		end
	end
	
	return push!(C, w)
end

__push_if_allowed(C::AbstractArray{T}, w::T, d::Integer) where T = __push_if_allowed!(copy(C), w, d)

#=
Convert inner-most elements into symbols
=#
__deepsym(a) = Symbol.(a)
__deepsym(a::AbstractArray) = __deepsym.(a)

__deepeltype(a) = __deepeltype(typeof(a))
(__deepeltype(::Type{T}) where T <: AbstractArray) = __deepeltype(eltype(T))
__deepeltype(::Type{T}) where T = T

__ensure_symbolic(Σ::AbstractArray) = __deepeltype(Σ) isa Symbol ? Σ : __deepsym(Σ)

__lessthanorequal(x, y)::Bool = isequal(x, y) || isless(x, y)

#=
Checks if a matrix has an identity in it.
e.g., [1 0 0 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0] does
=#
__has_identity(M::Matrix)::Bool = isequal(M[:, 1:size(M, 1)], I(size(M, 1))) ? true : false

#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
	
using LinearAlgebra: I

"""
	displaymatrix(M::AbstractArray)
	
Displays a matrix `M` in a compact form from the terminal.
"""
function displaymatrix(M::AbstractArray)
    return show(IOContext(stdout, :limit => true, :compact => true, :short => true), "text/plain", M); print("\n")
end

"""
	allequal_length(A::AbstractArray) -> Bool
	allequal_length(a, b...) -> Bool

Check that all elements in a list are of equal length.
"""
@inline function allequal_length(A::AbstractArray)
    length(A) < 2 && return true
	
    @inbounds for i in 2:length(A)
        isequal(length(A[1]), length(A[i])) || return false
    end
	
    return true
end

@inline function allequal_length(a...)::Bool
    A = [a...]
    return allequal_length(A)
end

"""
	allequal(A::AbstractArray) -> Bool
	allequal(a, b...) -> Bool

Check that all elements in a list are equal to each other.
"""
@inline function allequal(A::AbstractArray)
    length(A) < 2 && return true
    
    @inbounds for i in 2:length(A)
        A[1] ≠ A[i] && return false
    end
    
    return true
end

@inline function allequal(a...)::Bool
    A = [a...]
    return allequal(A)
end

"""
	aredistinct(A::AbstractArray) -> Bool
	aredistinct(a, b...) -> Bool

Check that all elements in a list are distinct from every other element in the list.
"""
@inline function aredistinct(A::AbstractArray)
    length(A) < 2 && return true
    
    while ! iszero(length(A))
        a = pop!(A)
        a ∉ A || return false
    end
    
    return true
end

@inline function aredistinct(a...)::Bool
    A = [a...]
    return aredistinct(A)
end

"""
	arelessthan(x::Number, A::AbstractArray) -> Bool
	arelessthan(x::Number, a, b...) -> Bool

Check that all elements in a list are less than a given x.
"""
@inline function arelessthan(x::Number, A::AbstractArray)
    @inbounds for a in A
        a < x || return false
    end
    
    return true
end

@inline function arelessthan(x::Number, a::Number...)::Bool
    A = [a...]
    return arelessthan(x, A)
end

"""
	push_if_allowed!(C::AbstractArray{T}, w::T, d::Integer)

Takes in an array and a word.  As long as the word does not mean that the distance is smaller than d, we add w to the array.  *This is a mutating function.  Use `push_if_allowed` for a non-mutating version of this function.*
"""
function push_if_allowed!(C::AbstractArray{T}, w::T, d::Integer) where T
	isempty(C) && return push!(C, w)
	
	for c in C
		if hamming_distance(c, w) < d
			return nothing
		end
	end
	
	return push!(C, w)
end

push_if_allowed(C::AbstractArray{T}, w::T, d::Integer) where T = push_if_allowed!(copy(C), w, d)

"""
	deepsym(a::AbstractArray)

Convert inner-most elements into symbols
"""
deepsym(a) = Symbol.(a)
deepsym(a::AbstractArray) = deepsym.(a)

"""
	deepeltype(a::AbstractArray) -> Type

Returns the type of the inner-most element in a nested array structure.
"""
deepeltype(a) = deepeltype(typeof(a))
(deepeltype(::Type{T}) where T <: AbstractArray) = deepeltype(eltype(T))
deepeltype(::Type{T}) where T = T

"""
	ensure_symbolic!(Σ::AbstractArray) -> AbstractArray

Ensures that the inner-most elements of a nested array structure are of the type `Symbol`.  *This is a mutating function.  Use its twin, non-mutating function, `ensure_symbolic`, if you need a non-mutating version of this.*
"""
ensure_symbolic!(Σ::AbstractArray) = deepeltype(Σ) isa Symbol ? Σ : deepsym(Σ)

"""
	ensure_symbolic(Σ::AbstractArray) -> AbstractArray

Ensures that the inner-most elements of a nested array structure are of the type `Symbol`.
"""
ensure_symbolic(Σ::AbstractArray) = ensure_symbolic(copy(Σ))

"""
	has_identity(M::Matrix) -> Bool

Checks if a matrix has an identity in it.
	
Examples:

	julia> A = [1 0 0 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]
	3×6 Array{Int64,2}:
	 1  0  0  2  3  0
	 0  1  0  1  2  2
	 0  0  1  4  3  0

	julia> B = [1 0 1 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]
	3×6 Array{Int64,2}:
	 1  0  1  2  3  0
	 0  1  0  1  2  2
	 0  0  1  4  3  0

	julia> has_identity(A)
	true

	julia> has_identity(B)
	false
"""
has_identity(M::Matrix) = isequal(M[:, 1:size(M, 1)], I(size(M, 1))) ? true : false

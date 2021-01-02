#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

"""
	displaymatrix(M::AbstractArray)
	
Displays a matrix `M` in a compact form from the terminal.
"""
function displaymatrix(M::AbstractArray)
    return show(IOContext(stdout, :limit => true, :compact => true, :short => true), "text/plain", M); print("\n")
end

_Iterable{T} = Union{AbstractArray{T}, NTuple{N, T}} where N

"""
	allequal_length(A) -> Bool
	allequal_length(a, b...) -> Bool

Check that all elements in a list are of equal length.
"""
@inline function allequal_length(A::_Iterable{T}) where T
    length(A) < 2 && return true
	
    @inbounds for i in 2:length(A)
        isequal(length(first(A)), length(A[i])) || return false
    end
	
    return true
end
@inline allequal_length(a::T...) where {T} = allequal_length(a)

"""
	allequal(A) -> Bool
	allequal(a, b...) -> Bool

Check that all elements in a list are equal to each other.
"""
@inline function allequal(A::_Iterable{T}) where T
    length(A) < 2 && return true
    
    @inbounds for i in 2:length(A)
        first(A) ≠ A[i] && return false
    end
    
    return true
end
@inline allequal(a::T...) where {T} = allequal(a)

"""
	aredistinct(A) -> Bool
	aredistinct(a, b...) -> Bool

Check that all elements in a list are distinct from every other element in the list.
"""
@inline function aredistinct(A::_Iterable{T}) where T
    length(A) < 2 && return true
    
    while ! iszero(length(A))
        a = pop!(A)
        a ∉ A || return false
    end
    
    return true
end
@inline aredistinct(a::T...) where {T} = aredistinct(a)

"""
	arelessthan(x::Number, A) -> Bool
	arelessthan(x::Number, a, b...) -> Bool

Check that all elements in a list are less than a given x.
"""
@inline function arelessthan(x::Number, A::_Iterable{T}) where T
    @inbounds for a in A
        a < x || return false
    end
    
    return true
end
@inline arelessthan(x::Number, a::Number...) = arelessthan(x, a)

"""
	areequalto(x::Number, A::AbstractArray) -> Bool
	areequalto(x::Number, a, b...) -> Bool

Check that all elements in a list are equal to a given x.
"""
@inline function areequalto(x, A::_Iterable{T}) where T
    @inbounds for a in A
        a != x || return false
    end
    
    return true
end
@inline areequalto(x, a::Number...) = areequalto(x, a)

"""
	deepsym(a::AbstractArray)

Convert inner-most elements into symbols
"""
deepsym(a) = Symbol.(a)
deepsym(a::_Iterable{T}) where {T} = deepsym.(a)

"""
	deepeltype(a::AbstractArray) -> Type

Returns the type of the inner-most element in a nested array structure.
"""
deepeltype(a) = deepeltype(typeof(a))
deepeltype(::Type{T}) where {T <: _Iterable{R}} where {R} = deepeltype(eltype(T))
deepeltype(::Type{T}) where T = T

"""
	ensure_symbolic!(Σ) -> typeof(Σ)

Ensures that the inner-most elements of a nested array structure are of the type `Symbol`.  *This is a mutating function.  Use its twin, non-mutating function, `ensure_symbolic`, if you need a non-mutating version of this.*
"""
ensure_symbolic!(Σ) = deepeltype(Σ) isa Symbol ? Σ : deepsym(Σ)

"""
	ensure_symbolic(Σ) -> typeof(Σ)

Ensures that the inner-most elements of a nested array structure are of the type `Symbol`.
"""
ensure_symbolic(Σ) = ensure_symbolic!(copy(Σ))

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
has_identity(M::Matrix) = isequal(M[:, axes(M, 1)], I(size(M, 1))) ? true : false

"""
	sizeof_perfect_code(q::Number, n::Number, d::Number) -> Number

Calculates the number of gigabytes required to store a perfect code of parameters q, n, and d.
"""
function sizeof_perfect_code(q::Int, n::Int, d::Int)
	return (sizeof(ntuple(_ -> gensym(), n)) * hamming_bound(big.([q, n, d])...)) / (2^30)
end
sizeof_perfect_code(q::Number, n::Number, d::Number) = sizeof_perfect_code(round.(BigInt, [q, n, d])...)

"""
	sizeof_perfect_code(q::Number, n::Number) -> Number

Calculates the number of gigabytes required to store all unique words of length n from an alphabet of size q.
"""
function sizeof_all_words(q::Int, n::Int)
	return (sizeof(ntuple(_ -> gensym(), n)) * big(q)^n) / (2^30)
end
sizeof_all_words(q::Number, n::Number) = sizeof_all_words(round.(BigInt, (q, n))...)

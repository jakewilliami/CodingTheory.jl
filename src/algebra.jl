#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "utils.jl"))
include(joinpath(dirname(@__FILE__), "rref.jl"))

using LinearAlgebra: I

"""
	normal_form!(M::AbstractArray, n::Integer) -> Matrix{Integer}

Convert a matrix M into normal form under modulo n via Gauss-Jordan elimination.  *This directly changes the matrix M.  Use `normal_form` for a non-mutating version of this function.*

Parameters:
  - M::AbstractArray: A matrix of integers.
  - n::Integer: The modulus of the finite field.
  
Returns:
  - Matrix{Integer}: A matrix in normal form from Gauss-Jordan elimination.
"""
function normal_form!(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	return rref!(M, n, colswap=false)
end

"""
	normal_form(M::AbstractArray, n::Integer) -> Matrix{Integer}

Convert a matrix M into normal form under modulo n via Gauss-Jordan elimination.

Parameters:
  - M::AbstractArray: A matrix of integers.
  - n::Integer: The modulus of the finite field.
  
Returns:
  - Matrix{Integer}: A matrix in normal form from Gauss-Jordan elimination.
"""
function normal_form(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	normal_form!(copy(M), n)
end

"""
	equivalent_code!(M::AbstractArray, n::Integer) -> Matrix{Integer}

Peforms Gauss-Jordan elimination on a matrix M, but allows for column swapping.  This constructs an "equivalent" matrix.  *This directly changes the matrix M.  Use `equivalent_code` for a non-mutating version of this function.*

Parameters:
  - M::AbstractArray: A matrix of integers.
  - n::Integer: The modulus of the finite field.
  
Returns:
  - Matrix{Integer}: A which represents an "equivalent" code to that of the matrix M.
"""
function equivalent_code!(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	return rref!(M, n, colswap=true)
end

"""
	equivalent_code(M::AbstractArray, n::Integer) -> Matrix{Integer}

Peforms Gauss-Jordan elimination on a matrix M, but allows for column swapping.

Parameters:
  - M::AbstractArray: A matrix of integers.
  - n::Integer: The modulus of the finite field.
  
Returns:
  - Matrix{Integer}: A which represents an "equivalent" code to that of the matrix M.
"""
function equivalent_code(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	equivalent_code!(copy(M), n)
end

"""
	generator!(M::AbstractArray, n::Integer; colswap::Bool=true) -> Matrix{Integer}

Constructs a generator matrix of the code, depending on if you allow for column swapping or not.   This function uses `normal_form!` or `equivalent_code!`.  *This directly changes the matrix M.  Use `generator` for a non-mutating version of this function.*

Parameters:
  - M::AbstractArray: A matrix of integers.
  - n::Integer: The modulus of the finite field.
  - colswap::Bool (kwarg): A boolean flag indicating whether or not you allow for swapping of columns when constructing the generating matrix.
	
Returns:
  - Matrix{Integer}: A generating matrix.
"""
function generator!(M::AbstractArray{T}, n::Integer; colswap::Bool=false)::Matrix{T} where T <: Integer
	return colswap ? equivalent_code!(M, n) : normal_form!(M, n)
end

"""
	generator(M::AbstractArray, n::Integer; colswap::Bool=true) -> Matrix{Integer}

Constructs a generator matrix of the code, depending on if you allow for column swapping or not.  This function uses `normal_form` or `equivalent_code`.

Parameters:
  - M::AbstractArray: A matrix of integers.
  - n::Integer: The modulus of the finite field.
  - colswap::Bool (kwarg): A boolean flag indicating whether or not you allow for swapping of columns when constructing the generating matrix.
	
Returns:
  - Matrix{Integer}: A generating matrix.
"""
function generator(M::AbstractArray{T}, n::Integer; colswap::Bool=true)::Matrix{T} where T <: Integer
	generator!(copy(M), n, colswap=colswap)
end

"""
	parity_check(M::AbstractArray, n::Integer) -> Matrix{Integer}
	
Constructs a parity check matrix.  This is calculated from taking the non-identity part of a matrix in normal form (or equivalent &mdash; see `generator`), transposing it, multiplying it by negative one, and appending to it an appropriate sized identity matrix.
	
Parameters:
  - M::AbstractArray: A matrix of integers.
  - n::Integer: The modulus of the finite field.
	
Returns:
  - Matrix{Integer}: A parity check matrix.
"""
function parity_check(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	if ! has_identity(M)
		throw(error("This matrix is not in normal form.  Use normal_form or equivalent_code."))
	end
	
	minus_Dᵀ = mod.(-one(Integer) .* transpose(M[:, size(M, 1)+1:end]), n)
	return H = [minus_Dᵀ I(size(minus_Dᵀ, 1))]
end

"""
	syndrome(v̲::Vector, Hᵀ::AbstractArray, n::Integer) -> Matrix{Integer}
	
Calculates the syndrome of a given word v̲ and a parity check matrix, transposed (Hᵀ), under modulo n.

Parameters:
  - v̲::Vector: A word in the code.
  - Hᵀ::AbstractArray: The transpose of a parity check matrix.
  
Returns:
  - Vector: The syndrome of a word in the code.
  
Examples:
	
	julia> syndrome([0, 2, 1, 2, 0, 1, 0], [1 1 1; 1 0 2; 2 0 1; 1 1 2; 1 0 0; 0 1 0; 0 0 1], 3)
	1×3 Array{Int64,2}:
	0  0  0
	
	julia> syndrome([0, 2, 1, 2, 0, 1, 0], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3)
	1×3 Array{Int64,2}:
	 0  0  0
"""
function syndrome(v̲::Vector, Hᵀ::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	return mod.(v̲' * Hᵀ, n)
end

"""
	isincode(v̲::Vector, Hᵀ::AbstractArray{T}, n::Integer) -> Bool

If the syndrome of a code is the zero vector, then the word used to calculate the syndrome is in the code.

Parameters:
  - v̲::Vector: A word.
  - Hᵀ::AbstractArray: The transpose of a parity check matrix.
  
Returns:
  - Bool: If the word is in the code or not (true or false).
"""
function isincode(v̲::Vector, Hᵀ::AbstractArray{T}, n::Integer)::Bool where T <: Integer
	iszero(syndrome(v̲, Hᵀ, n))
end

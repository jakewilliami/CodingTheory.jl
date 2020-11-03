#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
	
include(joinpath(dirname(@__FILE__), "abstract_types.jl"))
include(joinpath(dirname(@__FILE__), "utils.jl"))

using Polynomials
using Mods
using LinearAlgebra

"""
	struct FinitePolynomial <: FiniteField
		
Has parameters `p`, which is an abstract polynomial, and `n` which is the modulus of the field under which the molynomial is defined.
"""
struct FinitePolynomial <: FiniteField
	p::AbstractPolynomial
	n::Integer
end

"""
	FinitePolynomial(p::AbstractPolynomial, n::Integer)

A constructor method for `FinitePolynomial`.  Takes in a polynomial `p` and a number `n`, and constructs a polynomial under modulo n.
"""
function FinitePolynomial(p::AbstractPolynomial, n::Integer)
	p = Polynomial(mod.(p.coeffs, n))
	new(p, n)
end

"""
	mod(p::Polynomial, n::Integer) -> Polynomial
	
Uses the FinitePolynomial constructor to return a polynomial `p` under modulus `n`.

Parameters:
  - p::Polynomial: The input polynomial.
  - n::Integer: The modulus of the field.
  
Returns
  - Polynomial: A polynomial modulo n.
"""
Base.mod(p::Polynomial, n::Integer) = FinitePolynomial(p, n).p

"""
	Polynomial(A::Union{Tuple, AbstractArray}, n::Integer) -> Polynomial
	
Constructs a polynomial under modulo `n`.

Parameters:
  - A::Union{Tuple, AbstractArrau}: The polynomial coefficients.
  - n::Integer: The modulus of the field.
  
Returns
  - Polynomial: A polynomial modulo n.
"""
Polynomial(A::Union{Tuple, AbstractArray}, n::Integer) = mod(Polynomial(A), n)

"""
	list_polys(n::Integer, m::Integer) -> Array
	
Lists all polynomials of degree less than to `n` under modulo `m`.

Parameters:
  - n::Integer: Highest degree of polynomial.
  - m::Integer: The modulus of the field.
  
Returns:
  - Array: An array of polynomials of degree less than n, under modulo m.
"""
function list_polys(n::Integer, m::Integer)
	return collect(Polynomial(collect(t)) for t in Iterators.product([0:(m-1) for i in 1:n]...))
end

"""
	multiplication_table(degree::Integer, modulo::Integer) -> Matrix
	
Returns a table (matrix) of the multiplication of all combinations of polynomials for degree less than `degree`, under modulo `modulo`.

Parameters:
  - degree::Integer: Highest degree of polynomial.
  - modulo::Integer: The modulus of the field.
  
Returns:
  - Matrix: A multiplication table of all polynomials with degree less than n, under modulus.
"""
function multiplication_table(degree::Integer, modulo::Integer)
	polys = list_polys(degree, modulo)
	number_of_polys = length(polys)
	poly_matrix = Matrix{Polynomial}(undef, number_of_polys, number_of_polys)
	
	for i in 1:number_of_polys, j in 1:number_of_polys
		poly_matrix[i,j] = mod(polys[i]*polys[j], modulo)
	end

	return poly_matrix
end

"""
	list_span(u̲::Vector, v̲::Vector, modulo::Integer) -> Array

Given two vectors `u̲` and `v̲`, prints all linear combinations of those vectors, modulo `modulo`.  NB: this function can take in another vector, but is not yet generalised to more than three.

Parameters:
  - u̲::Vector: One vector.
  - v̲::Vector: Another vector.
  - modulo::Integer: The modulus of the field.
  
Returns:
  - Array: All vectors in the span of u̲ and v̲, under modulo.
"""
function list_span(u̲::Vector, v̲::Vector, modulo::Integer)
	span = Vector[]
	
	for λ in 0:modulo-1, γ in 0:modulo-1
		w̲ = mod.(λ*u̲ + γ*v̲, modulo)
		if w̲ ∉ span
			push!(span, w̲)
		end
	end
	
	return span
end

function list_span(u̲::Vector, v̲::Vector, t̲::Vector, modulo::Integer)
	span = Vector[]
	
	for λ in 0:modulo-1, γ in 0:modulo-1, α in 0:modulo-1
		w̲ = mod.(λ*u̲ + γ*v̲ + α*t̲, modulo)
		if w̲ ∉ span
			push!(span, w̲)
		end
	end
	
	return span
end

"""
	islinear(C::Vector, modulo::Integer; verbose::Bool=false) -> Bool
	
Determines whether a code `C` is a linear code (i.e., if it is closed under addition, scalar multiplication, and has the zero vector in it).

Parameters:
  - C::Vector: A code, typically consisting of multiple vectors or strings.
  - modulo::Integer: The modulus of the field under which you are working.
  - verbose::Bool (kwarg): print the point at which C fails, if it does.
  
Returns:
  - Bool: Whether or not the code `C` is linear (true or false).
"""
function islinear(C::Vector, modulo::Integer; verbose::Bool=false)
	allequal_length(C) || return false # not all codes are of the same length
	block_length = length(C[1])
	𝟎 = fill(0, block_length)
		
	if 𝟎 ∉ C
		if verbose
			println("The zero vector 0̲ is not in C.\n")
		end
		return false # the zero vector is not in the code
	end
	
	for c̲ ∈ C
		for λ in 0:modulo-1
			if mod.(λ*c̲, modulo) ∉ C
				if verbose
					println(λ, " ⋅ ", c̲, " = ", mod.(λ*c̲, modulo), " ∉ C\n")
				end
				return false # this code isn't closed under scalar multiplication
			end
		end
		
		for c̲ ∈ C, c̲′ ∈ C
			if c̲ ≠ c̲′
				if mod.(c̲ + c̲′, modulo) ∉ C
					if verbose
						println(c̲, " + ", c̲′, " = ", mod.(c̲ + c̲′, modulo), " ∉ C\n")
					end
					return false # this code isn't closed under addition
				end
			end
		end
	end
	
	if verbose
		println()
	end
	
	return true
end

"""
	isirreducible(f::AbstractPolynomial, modulo::Integer) -> Bool
	
Checks if a polynomial is irreducible.
	
Parameters:
  - f::Polynomial: The polynomial you need to check.
  - modulo::Integer: The modulus under which you are working.
  
Returns:
  - Bool: Whether or not the polynomial is irreducible (true or false).
"""
function isirreducible(f::AbstractPolynomial, modulo::Integer)
	deg = length(f.coeffs) - 1
	f = mod(f, deg)
	polys = list_polys(deg, modulo)

	for a in polys, b in polys
		isequal(f, mod(a*b, modulo)) && return false
	end
	
	return true
end

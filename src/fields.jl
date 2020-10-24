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

struct FinitePolynomial <: FiniteField
	p::AbstractPolynomial
	n::Integer
	
	function FinitePolynomial(p::AbstractPolynomial, n::Integer)
		p = Polynomial(mod.(p.coeffs, n))
		new(p, n)
	end
end

Base.mod(p::Polynomial, n::Integer) = FinitePolynomial(p, n).p
Polynomial(A::Union{Tuple, AbstractArray}, n::Integer) = mod(Polynomial(A), n)

function list_polys(n::Integer, m::Integer)::AbstractArray
	return collect(Polynomial(collect(t)) for t in Iterators.product([0:(m-1) for i in 1:n]...))
end

function multiplication_table(degree::Integer, modulo::Integer)::Matrix
	polys = list_polys(degree, modulo)
	number_of_polys = length(polys)
	poly_matrix = Matrix{Polynomial}(undef, number_of_polys, number_of_polys)
	
	for i in 1:number_of_polys, j in 1:number_of_polys
		poly_matrix[i,j] = mod(polys[i]*polys[j], modulo)
	end

	return poly_matrix
end

function list_span(u̲::Vector, v̲::Vector, modulo::Integer)::Array{Array{Int, 1}}
	span = Vector[]
	
	for λ in 0:modulo-1, γ in 0:modulo-1
		w̲ = mod.(λ*u̲ + γ*v̲, modulo)
		if w̲ ∉ span
			push!(span, w̲)
		end
	end
	
	return span
end

function list_span(u̲::Vector, v̲::Vector, t̲::Vector, modulo::Integer)::Array{Array{Int, 1}}
	span = Vector[]
	
	for λ in 0:modulo-1, γ in 0:modulo-1, α in 0:modulo-1
		w̲ = mod.(λ*u̲ + γ*v̲ + α*t̲, modulo)
		if w̲ ∉ span
			push!(span, w̲)
		end
	end
	
	return span
end

function islinear(C::Vector, modulo::Integer; verbose::Bool=false)::Bool
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

function isirreducible(f::AbstractPolynomial, modulo::Integer)::Bool
	deg = length(f.coeffs) - 1
	f = mod(f, deg)
	polys = list_polys(deg, modulo)

	for a in polys, b in polys
		isequal(f, mod(a*b, modulo)) && return false
	end
	
	return true
end

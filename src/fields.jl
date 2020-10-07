#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

using Polynomials
using Mods
using LinearAlgebra
using RowEchelon

Base.mod(p::Polynomial, n::Integer) = Polynomial(mod.(p.coeffs, n))
RowEchelon.rref(M::AbstractArray, n::Integer) = mod.(Int.(RowEchelon.rref(M)), n)
RowEchelon.rref!(M::AbstractArray, n::Integer) = mod.(Int.(RowEchelon.rref!(M)), n)

@generated function list_polys(::Val{n}, m::Integer)::AbstractArray where {n}
    quote
		polys = Polynomial[]
		
        Base.Cartesian.@nloops $n i d -> 0:m-1 begin
			push!(polys, Polynomial([(Base.Cartesian.@ntuple $n i)...]))
        end
		
		return polys
    end
end

list_polys(n::Integer, modulo::Integer)::AbstractArray = list_polys(Val(n), modulo)

function multiplication_table(degree::Integer, modulo::Integer)::Matrix
	polys = list_polys(degree, modulo)
	number_of_polys = length(polys)
	poly_matrix = Matrix{Polynomial}(undef, number_of_polys, number_of_polys)
	
	for i in 1:number_of_polys, j in 1:number_of_polys
		poly_matrix[i,j] = mod(polys[i]*polys[j], modulo)
	end

	return poly_matrix
end

function list_span(u̲::Vector, v̲::Vector, modulo::Integer)::Vector
	span = Vector[]
	
	for λ in 0:modulo-1, γ in 0:modulo-1
		w̲ = mod.(λ*u̲ + γ*v̲, modulo)
		if w̲ ∉ span
			push!(span, w̲)
		end
	end
	
	return span
end

function list_span(u̲::Vector, v̲::Vector, t̲::Vector, modulo::Integer)::Vector
	span = Vector[]
	
	for λ in 0:modulo-1, γ in 0:modulo-1, α in 0:modulo-1
		w̲ = mod.(λ*u̲ + γ*v̲ + α*t̲, modulo)
		if w̲ ∉ span
			push!(span, w̲)
		end
	end
	
	return span
end

@inline function _allequal_length_(A::Array)::Bool
    length(A) < 2 && return true
	
    e1 = A[1]
	
    @inbounds for i in 2:length(A)
        isequal(length(A[i]), length(e1)) || return false
    end
	
    return true
end

function islinear(C::Vector, modulo::Integer)#::Bool
	_allequal_length_(C) || return false # not all codes are of the same length
	block_length = length(C[1])
	𝟎 = fill(0, block_length)
		
	𝟎 ∈ C || return false # the zero vector is not in the code
	
	for c̲ ∈ C
		for λ in 0:modulo-1
			mod.(λ*c̲, modulo) ∈ C || return false # this code isn't closed under scalar multiplication
		end
		
		for c̲ ∈ C, c̲′ ∈ C
			if c̲ ≠ c̲′
				mod.(c̲ + c̲′, modulo) ∈ C || return false # this code isn't closed under addition
			end
		end
	end
	
	return true
end

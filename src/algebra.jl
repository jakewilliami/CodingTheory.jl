#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "utils.jl"))
include(joinpath(dirname(@__FILE__), "rref.jl"))

using LinearAlgebra: I

function normal_form!(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	return rref!(M, n, colswap=false)
end

(normal_form(M::AbstractArray{T},
	n::Integer
	)::Matrix{T}) where T <: Integer = normal_form!(copy(M), n)
	
function equivalent_code!(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	return rref!(M, n, colswap=true)
end

(equivalent_code(M::AbstractArray{T},
	n::Integer
	)::Matrix{T}) where T <: Integer = equivalent_code!(copy(M), n)
	
function generator!(M::AbstractArray{T}, n::Integer; colswap::Bool=false)::Matrix{T} where T <: Integer
	return colswap ? equivalent_code!(M, n, colswap=colswap) : normal_form!(M, n, colswap=colswap)
end

(generator(M::AbstractArray{T},
	n::Integer;
	colswap::Bool=true
	)::Matrix{T}) where T <: Integer = generator!(copy(M), n, colswap=colswap)
	
function parity_check(M::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	if ! __has_identity(M)
		throw(error("This matrix is not in normal form.  Use normal_form or equivalent_code."))
	end
	
	minus_Dᵀ = mod.(-one(Integer) .* transpose(M[:, size(M, 1)+1:end]), n)
	return H = [minus_Dᵀ I(size(minus_Dᵀ, 1))]
end

function syndrome(v̲::Vector, Hᵀ::AbstractArray{T}, n::Integer)::Matrix{T} where T <: Integer
	return mod.(v̲' * Hᵀ, n)
end

(isincode(v̲::Vector,
	Hᵀ::AbstractArray{T},
	n::Integer
	)::Bool) where T <: Integer = iszero(syndrome(v̲, Hᵀ, n))

  
#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

#=
Convert inner-most elements into symbols
=#
__deepsym(a) = Symbol.(a)
__deepsym(a::AbstractArray) = __deepsym.(a)

__lessthanorequal(x, y)::Bool = isequal(x, y) || isless(x, y)

function hamming_distance(w₁, w₂)::Integer
    if ! isequal(length(w₁), length(w₂))
        throw(error("Cannot compute Hamming Distance on strings of unequal length."))
    end
    
    distance = 0
        
    for i in 1:length(w₁)
        if w₁[i] ≠ w₂[i]
            distance += 1
        end
    end
    
    return distance
end

#=
Get the codewords of radius e of a ball centered at w
=#
function __hamming_space(relation::Function, Σⁿ::AbstractArray{T}, w::AbstractArray, e::Integer) where T# <: AbstractArray
	e < 0 && throw(error("e (the ball \"radius\") must be a non-negative number."))
	
	if eltype(w) isa Symbol
	else
		w = __deepsym(w)
	end
	
	ball = []
	
	for v in __deepsym(Σⁿ)
		relation(hamming_distance(w, v), e) && push!(ball, v)
	end
	
	return ball
end

hamming_ball(Σⁿ::AbstractArray{T}, w::AbstractArray, e::Integer) where T =
	__hamming_space(__lessthanorequal, Σⁿ, w, e)
hamming_sphere(Σⁿ::AbstractArray{T}, w::AbstractArray, e::Integer) where T =
	__hamming_space(isequal, Σⁿ, w, e)

#=
Convert contents of nested arrays into symbols
=#
function code_distance(C::AbstractArray{T})::Integer where T <: AbstractArray{Any}
	return code_distance(__deepsym(C))
end

# function code_distance(C::AbstractArray{Int})::Integer
# 	return code_distance(collect(copy.(eachcol(A))))
# end

#=
Find the minimum hamming distance in the code for all unique letters
=#
function code_distance(C::AbstractArray{T})::Integer where T
	distances = []
	
	for c in C, c′ in C
		if c ≠ c′
			push!(distances, hamming_distance(c, c′))
		end
	end

	return minimum(distances)
end

#=
A wrapper to get the code distance after pushing a word to the code.
=#
function code_distance!(C::AbstractArray{T}, w::T)::Integer where T
	push!(C, w)
	
	return code_distance(C)
end

(code_distance(C::AbstractArray{T}, w::T)::Integer) where T = code_distance!(copy(C), w)

function t_error_detecting(C::AbstractArray{T}, t::Integer)::Bool where T <: AbstractArray{Int}
	code_distance(C) ≥ t + 1 && return true
	return false
end

function t_error_correcting(C::AbstractArray{T}, t::Integer)::Bool where T <: AbstractArray{Int}
	code_distance(C) ≥ 2*t + 1 && return true
	return false
end

function find_error_detection_max(C::AbstractArray{T}, modulo::Integer)::Integer where T <: AbstractArray{Int}
	for t in modulo-1:-1:0
		t_error_detecting(C, t) && return t
	end
end

function find_error_correction_max(C::AbstractArray{T}, modulo::Integer)::Integer where T <: AbstractArray{Int}
	for t in modulo-1:-1:0
		t_error_correcting(C, t) && return t
	end
end

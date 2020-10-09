  
#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

function hamming_distance(s1::Union{AbstractString, AbstractArray}, s2::Union{AbstractString, AbstractArray})::Integer
    if ! isequal(length(s1), length(s2))
        throw(error("Cannot compute Hamming Distance on strings of unequal length."))
    end
    
    distance = 0
        
    for i in 1:length(s1)
        if s1[i] ≠ s2[i]
            distance += 1
        end
    end
    
    return distance
end

function hamming_ball(Σⁿ::AbstractArray{T}, w::Array{Int, 1}, e::Integer)::Array{Array{Int, 1}} where T <: AbstractArray{Int}
	e > 0 || throw(error("e (the ball \"radius\") must be a non-negative number."))
	
	ball = Vector[]
	
	for v in Σⁿ
		hamming_distance(w, v) ≤ e && push!(ball, v)
	end
	
	return ball
end

function code_distance(C::AbstractArray{String})::Integer
	distances = Integer[]
	
	for c in C, c′ in C
		if c ≠ c′
			push!(distances, hamming_distance(c, c′))
		end
	end

	return minimum(distances)
end

function code_distance(C::AbstractArray{T})::Integer where T <: AbstractArray{Int}
	string_code = String[]
	
	for c in C
		push!(string_code, join(c))
	end
	
	return code_distance(string_code)
end

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

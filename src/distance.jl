  
#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

"""
	hamming_distance(w₁, w₂) -> Integer

The Hamming distance of two words is the number of changes that need to be made to each letter in the word for the words to be the same.  This does not work for words of unequal length.
	
Parameters:
  - w₁: A word.
  - w₂: Another word.
  
Returns:
  - Integer: the number of changes needing to be made to one word for it to be identical to the other.
"""
function hamming_distance(w₁, w₂)
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

function __hamming_space(relation::Function, Σⁿ::AbstractArray{T}, w::AbstractArray, e::Integer) where T# <: AbstractArray
	e < 0 && throw(error("e (the ball/sphere \"radius\") must be a non-negative number."))
	
	w = ensure_symbolic(w)
	ball = []
	
	for v in deepsym(Σⁿ)
		relation(hamming_distance(w, v), e) && push!(ball, v)
	end
	
	return ball
end

"""
	hamming_ball(Σⁿ::AbstractArray, w::AbstractArray, e::Integer) -> Vector{Vector}

Get the codewords of radius e of a ball centered at word w.  That is, all words whose distance from w is less than or equal to the radius.

Parameters:
  - Σⁿ::AbstractArray: An array of words in the code.
  - w::AbstractArray: A word.
  - e::Integer: The radius of the ball.
  
Returns:
  - AbstractArray: The list of words in Σⁿ whose distance from w is less than or equal to e.  Returns an array of array of symbols.
"""
hamming_ball(Σⁿ::AbstractArray{T}, w::AbstractArray, e::Integer) where T =
	__hamming_space(≤, Σⁿ, w, e)

"""
	hamming_sphere(Σⁿ::AbstractArray, w::AbstractArray, e::Integer) -> Vector{Vector}

Get the codewords of radius e of a sohere centered at word w.  That is, all words whose distance from w is exactly equal to to the radius.

Parameters:
  - Σⁿ::AbstractArray: An array of words in the code.
  - w::AbstractArray: A word.
  - e::Integer: The radius of the ball.
  
Returns:
  - AbstractArray: The list of words in Σⁿ whose distance from w is exactly equal to e.  Returns an array of array of symbols.
"""
hamming_sphere(Σⁿ::AbstractArray{T}, w::AbstractArray, e::Integer) where T =
	__hamming_space(isequal, Σⁿ, w, e)

"""
	code_distance(C::AbstractArray) -> Integer
	
Finds the distance of the code.  That is, given a code C, finds the minimum distance between any two words in the code, which are not the same. (Find the minimum hamming distance in the code for all unique letters).

Parameters:
  - C::AbstractArray: An array of words in the code.
  
Return:
  - Integer: the distance of the code.
"""
function code_distance(C::AbstractArray{T}) where T <: AbstractArray{Any}
	return code_distance(deepsym(C))
end

# function code_distance(C::AbstractArray{Int})::Integer
# 	return code_distance(collect(copy.(eachcol(A))))
# end

function code_distance(C::AbstractArray{T}) where T
	isempty(C) && return nothing
	distances = []
	
	for c in C, c′ in C
		if c ≠ c′
			push!(distances, hamming_distance(c, c′))
		end
	end

	return minimum(distances)
end

"""
	code_distance!(C::AbstractArray{T}, w::T) -> Integer

A wrapper to get the code distance after pushing a word to the code.  *This directly changes the matrix M.  Use `code_distance` for a non-mutating version of this function.*

Parameters:
  - C::AbstractArray: An array of words in the code.
  - w: A word to be appended to C.
  
Returns:
  - Integer: The code distance after adding w to C.
"""
function code_distance!(C::AbstractArray{T}, w::T) where T
	push!(C, w)
	return code_distance(C)
end

"""
	code_distance(C::AbstractArray{T}, w::T) -> Integer

A wrapper to get the code distance after pushing a word to the code.

Parameters:
  - C::AbstractArray: An array of words in the code.
  - w: A word to be appended to C.
  
Returns:
  - Integer: The code distance after adding w to C.
"""
function code_distance(C::AbstractArray{T}, w::T) where T
	code_distance!(copy(C), w)
end

"""
	t_error_detecting(C::AbstractArray{T}, t::Integer) -> Bool
	
Check if a given code C can detect t many errors.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - t::Integer: The number of errors you want to check that the code can detect.
  
Returns:
  - Bool: Yes, C can detect t errors, or no it cannot (true of false).
"""
function t_error_detecting(C::AbstractArray{T}, t::Integer) where T <: AbstractArray{Int}
	code_distance(C) ≥ t + 1 && return true
	return false
end

"""
	t_error_correcting(C::AbstractArray{T}, t::Integer) -> Bool
	
Check if a given code C can correct t many errors.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - t::Integer: The number of errors you want to check that the code can correct.
  
Returns:
  - Bool: Yes, C can correct t errors, or no it cannot (true of false).
"""
function t_error_correcting(C::AbstractArray{T}, t::Integer) where T <: AbstractArray{Int}
	code_distance(C) ≥ 2*t + 1 && return true
	return false
end


"""
	find_error_detection_max(C::AbstractArray{T}, modulo::Integer) -> Integer
	
Finds the greatest number t such that C is error t error detecting.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - moldulo::Integer: The modulus of the finite field.  The upper bound of t.
  
Returns:
  - Integer: The maximum number t such that the code is t error detecting.
"""
function find_error_detection_max(C::AbstractArray{T}, modulo::Integer) where T <: AbstractArray{Int}
	for t in modulo-1:-1:0
		t_error_detecting(C, t) && return t
	end
end

"""
	find_error_correction_max(C::AbstractArray{T}, modulo::Integer) -> Integer
	
Finds the greatest number t such that C is error t error correcting.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - moldulo::Integer: The modulus of the finite field.  The upper bound of t.
  
Returns:
  - Integer: The maximum number t such that the code is t error correcting.
"""
function find_error_correction_max(C::AbstractArray{T}, modulo::Integer) where T <: AbstractArray{Int}
	for t in modulo-1:-1:0
		t_error_correcting(C, t) && return t
	end
end

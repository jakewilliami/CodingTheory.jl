  
#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

"""
	hamming_distance(w₁, w₂) -> Int

The Hamming distance of two words is the number of changes that need to be made to each letter in the word for the words to be the same.  This does not work for words of unequal length.
	
Parameters:
  - w₁: A word.
  - w₂: Another word.
  
Returns:
  - Int: the number of changes needing to be made to one word for it to be identical to the other.
"""
function hamming_distance(w₁::T, w₂::T) where T <: AbstractWord
    if ! isequal(length(w₁), length(w₂))
        throw(error("Cannot compute Hamming Distance on strings of unequal length."))
    end
    
    distance = 0
    
	for (s₁, s₂) in zip(w₁, w₂)
		if s₁ ≠ s₂
			distance += 1
		end
	end
    
    return distance
end

function __hamming_space(relation::Function, Σⁿ::AbstractArray{T}, w::AbstractArray, e::Int) where T <: AbstractWord
	e < 0 && throw(error("e (the ball/sphere \"radius\") must be a non-negative number."))
	Σⁿ, w = deepsym(Σⁿ), ensure_symbolic(w)
	
	return eltype(Σⁿ)[v for v in deepsym(Σⁿ) if relation(hamming_distance(w, v), e)]
end

"""
	hamming_ball(Σⁿ::AbstractArray, w::Vector, e::Int) -> Vector{Vector}

Get the codewords of radius e of a ball centered at word w.  That is, all words whose distance from w is less than or equal to the radius.

Parameters:
  - Σⁿ::AbstractArray: An array of words in the code.
  - w::Vector: A word.
  - e::Int: The radius of the ball.
  
Returns:
  - AbstractArray: The list of words in Σⁿ whose distance from w is less than or equal to e.  Returns an array of array of symbols.
"""
hamming_ball(Σⁿ::AbstractArray{T}, w::Vector{S}, e::Int) where {T <: AbstractWord, S} =
	__hamming_space(≤, Σⁿ, w, e)

"""
	hamming_sphere(Σⁿ::AbstractArray, w::Vector, e::Int) -> Vector{Vector}

Get the codewords of radius e of a sohere centered at word w.  That is, all words whose distance from w is exactly equal to to the radius.

Parameters:
  - Σⁿ::AbstractArray: An array of words in the code.
  - w::Vector: A word.
  - e::Int: The radius of the ball.
  
Returns:
  - AbstractArray: The list of words in Σⁿ whose distance from w is exactly equal to e.  Returns an array of array of symbols.
"""
hamming_sphere(Σⁿ::AbstractArray{T}, w::Vector{S}, e::Int) where {T <: AbstractWord, S} =
	__hamming_space(isequal, Σⁿ, w, e)

"""
	code_distance(C::AbstractArray) -> Int
	
Finds the distance of the code.  That is, given a code C, finds the minimum distance between any two words in the code, which are not the same. (Find the minimum hamming distance in the code for all unique letters).

Parameters:
  - C::AbstractArray: An array of words in the code.
  
Return:
  - Int: the distance of the code.
"""
function code_distance(C::AbstractArray{T}) where T <: AbstractWord
	isempty(C) && return nothing
	C = deepsym(C)
	min_distance = nothing
	
	for c in C, c′ in C
		if c ≠ c′
			ham_dist = hamming_distance(c, c′)
			if isnothing(min_distance) || ham_dist < min_distance
				min_distance = ham_dist
			end
		end
	end
	
	return min_distance
end

"""
	code_distance!(C::AbstractArray{T}, w::T) -> Int

A wrapper to get the code distance after pushing a word to the code.  *This directly changes the matrix M.  Use `code_distance` for a non-mutating version of this function.*

Parameters:
  - C::AbstractArray: An array of words in the code.
  - w: A word to be appended to C.
  
Returns:
  - Int: The code distance after adding w to C.
"""
function code_distance!(C::AbstractArray{T}, w::T) where T <: AbstractWord
	push!(C, w)
	return code_distance(C)
end

"""
	code_distance(C::AbstractArray{T}, w::T) -> Int

A wrapper to get the code distance after pushing a word to the code.

Parameters:
  - C::AbstractArray: An array of words in the code.
  - w: A word to be appended to C.
  
Returns:
  - Int: The code distance after adding w to C.
"""
function code_distance(C::AbstractArray{T}, w::T) where T <: AbstractWord
	code_distance!(copy(C), w)
end

"""
	t_error_detecting(C::AbstractArray{T}, t::Int) -> Bool
	
Check if a given code C can detect t many errors.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - t::Int: The number of errors you want to check that the code can detect.
  
Returns:
  - Bool: Yes, C can detect t errors, or no it cannot (true of false).
"""
function t_error_detecting(C::AbstractArray{T}, t::Int) where T <: Union{AbstractArray{Int}, AbstractWord}
	code_distance(C) ≥ t + 1 && return true
	return false
end

"""
	t_error_correcting(C::AbstractArray{T}, t::Int) -> Bool
	
Check if a given code C can correct t many errors.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - t::Int: The number of errors you want to check that the code can correct.
  
Returns:
  - Bool: Yes, C can correct t errors, or no it cannot (true of false).
"""
function t_error_correcting(C::AbstractArray{T}, t::Int) where T <: Union{AbstractArray{Int}, AbstractWord}
	code_distance(C) ≥ 2*t + 1 && return true
	return false
end


"""
	find_error_detection_max(C::AbstractArray{T}, modulo::Int) -> Int
	
Finds the greatest number t such that C is error t error detecting.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - moldulo::Int: The modulus of the finite field.  The upper bound of t.
  
Returns:
  - Int: The maximum number t such that the code is t error detecting.
"""
function find_error_detection_max(C::AbstractArray{T}, modulo::Int) where T <: Union{AbstractArray{Int}, AbstractWord}
	for t in (modulo - 1):-1:0
		t_error_detecting(C, t) && return t
	end
end

"""
	find_error_correction_max(C::AbstractArray{T}, modulo::Int) -> Int
	
Finds the greatest number t such that C is error t error correcting.
	
Parameters:
  - C::AbstractArray: An array of words in the code.
  - moldulo::Int: The modulus of the finite field.  The upper bound of t.
  
Returns:
  - Int: The maximum number t such that the code is t error correcting.
"""
function find_error_correction_max(C::AbstractArray{T}, modulo::Int) where T <: Union{AbstractArray{Int}, AbstractWord}
	for t in modulo-1:-1:0
		t_error_correcting(C, t) && return t
	end
end

#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

include("distance.jl")
include("primes.jl")

"""
	rate(q::Integer, M::Integer, n::Integer) -> Real
	
Calculate the rate of a code.  That is, how efficient the code is.

Parameters:
  - q::Integer: the number of symbols in the code.
  - M::Integer: the size/number of elements in the code.
  - n::Integer: The word length.

Returns:
  - Real: Rate of the code.
"""
rate(q::Integer, M::Integer, n::Integer) = log(q, M) / n

__spheres(q::Integer, n::Integer, r::Integer) = sum([((big(q) - 1)^i) * binomial(big(n), big(i)) for i in 0:r])
__sphere_bound(round_func::Function, q::Integer, n::Integer, d::Integer) = round_func((big(q)^n) / __spheres(q, n, d))

"""
	sphere_covering_bound(q::Integer, n::Integer, d::Integer) -> Integer
	
Computes the sphere covering bound of a [n, d]q-code.

Parameters:
  - q::Integer: the number of symbols in the code.
  - n::Integer: the word length.
  - d::Integer: the distance of the code.
  
Returns:
  - Integer: the sphere covering bound.
"""
sphere_covering_bound(q::Integer, n::Integer, d::Integer) = __sphere_bound(ceil, q, n, d - 1)

"""
	sphere_packing_bound(q::Integer, n::Integer, d::Integer) -> Integer
	sphere_packing_bound(q::Integer, n::Integer, d::Integer, ::Rounding) -> Real
	
Computes the sphere packing bound of a [n, d]q-code.  The sphere packing bound is also known as the hamming bound.  You can use `hamming_bound` to compute the same thing.

Parameters:
  - q::Integer: the number of symbols in the code.
  - n::Integer: the word length.
  - d::Integer: the distance of the code.
  - ::Rounding: use the argument `no_round` in this position to preserve the rounding of the code &mdash; which usually by default rounds down.
  
Returns:
  - Integer: the sphere packing bound.
"""
sphere_packing_bound(q::Integer, n::Integer, d::Integer) = __sphere_bound(floor, q, n, Int(floor((d - 1) / 2)))
sphere_packing_bound(q::Integer, n::Integer, d::Integer, ::Rounding) = __sphere_bound(identity, q, n, Int(floor((d - 1) / 2)))
hamming_bound(q::Integer, n::Integer, d::Integer) = sphere_packing_bound(q, n, d)
hamming_bound(q::Integer, n::Integer, d::Integer, ::Rounding) = sphere_packing_bound(q, n, d, no_round)

"""
	sphere_packing_bound(q::Integer, n::Integer, d::Integer) -> Real
	
Computes the Singleton bound of a [n, d]q-code.

Parameters:
  - q::Integer: the number of symbols in the code.
  - n::Integer: the word length.
  - d::Integer: the distance of the code.
  
Returns:
  - Real: the Singleton bound.  Can round down, as it is an equivalent to the Hamming bound in that it is an upper bound.
"""
singleton_bound(q::Number, n::Number, d::Number) = float(big(q))^(big(n) - big(d) + 1)

"""
	construct_ham_matrix(r::Integer, q::Integer) -> Matrix
	
Construct a Hamming parity-check matrix.

Parameters:
  - r::Integer: number of rows of a parity check matrix.
  - q:::Integer: The size of the alphabet of the code.
  
Returns:
  - Matrix: The Hamming matrix, denoted as Ham(r, q)
"""
function construct_ham_matrix(r::Integer, q::Integer)::Matrix
    ncols = Int(floor((q^r - 1) / (q - 1)))
    M = Matrix{Integer}(undef, r, ncols)
    
    for i in 1:ncols
        M[:,i] = reverse(digits(parse(Int, string(i, base = q)), pad = r), dims = 1)
    end
    
    return M
end

"""
	isperfect(n::Integer, k::Integer, d::Integer, q::Integer) -> Bool
	
Checks if a code is perfect.  That is, checks if the number of words in the code is exactly the "Hamming bound", or the "Sphere Packing Bound".
	
Parameters:
  - q:::Integer: The size of the alphabet of the code.
  - n::Integer: The length of the words in the code (block length).
  - d::Integer: The distance of the code.
  - k::Integer: The dimension of the code.
  
Returns:
  - Bool: true or false
"""
function isperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)::Bool
	
	isprimepower(q) || throw(error("Cannot check if the code is perfect with q not a prime power."))
    
    M = q^k
    
    if isequal(sphere_packing_bound(q, n, d), M)
        return true
    else
        return false
    end
end

"""
	ishammingbound(r::Integer, q::Integer) -> Bool
	
Checks if the code is a perfect code that is of the form of a generalised Hamming code.
	
Parameters:
  - r::Integer: number of rows of a parity check matrix.
  - q:::Integer: The size of the alphabet of the code.
  
Returns:
  - Bool: true or false
"""
function ishammingperfect(r::Integer, q::Integer)::Bool
    n = 2^r - 1
    k = n - r
    M = q^k
    d = size(construct_ham_matrix(r, q))[1] # the number of rows of the hamming matrix (which is, by design, linearly independent)
    d = r
    # r is dim of dueal code; dim of code itself is block length minus r
    println(n)
    println((q^r - 1) / (q - 1))
    
    if isequal(n, (q^r - 1) / (q - 1)) && isequal(d, 3) && isequal(M, q^(((q^r - 1) / (q - 1)) - r))
        return true
    end
    
    return false
end

"""
	ishammingperfect(n::Integer, k::Integer, d::Integer, q::Integer) -> Bool
	ishammingperfect(q::Integer, n::Integer, d::Integer) -> Bool
	
Checks if the code is a perfect code that is of the form of a generalised Hamming code.
	
Parameters:
  - q:::Integer: The size of the alphabet of the code.
  - n::Integer: The length of the words in the code (block length).
  - d::Integer: The distance of the code.
  - k::Integer: The dimension of the code.
  
Returns:
  - Bool: true or false
"""
function ishammingperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)
    
    ! isprimepower(q) && return false
    
    M = q^k
    r = log(ℯ, ((n * log(ℯ, 1)) / (log(ℯ, 2))) + 1) / log(ℯ, 2)
        
    if isequal(n, (q^(r - 1)) / (q - 1)) && isequal(d, 3) && isequal(M, q^(((q^r - 1) / (q - 1)) - r))
        return true
    end
    
    return false
end

function ishammingperfect(q::Integer, n::Integer, d::Integer)
	! isprimepower(q) && return false # we are working in finite fields, so q must be a prime power
	d ≠ 3 && return false
	
	r = 1
	while ((q^r - 1) / (q - 1)) < n
		r = r + 1
	end
	
	return isequal(((q^r - 1) / (q - 1)), n) ? true : false
end

"""
	isgolayperfect(n::Integer, k::Integer, d::Integer, q::Integer) -> Bool
	
Golay found two perfect codes.  `isgolayperfect` checks if a code of block length n, distance d, alphabet size q, and dimension k, is a perfect code as described by Golay.

Parameters:
  - n::Integer: The block length of words in the code (e.g., word "abc" has block length 3).
  - k::Integer: The dimension of the code.
  - d::Integer: The distance of the code (i.e., the minimum distance between codewords in the code).
  - q::Integer: An integer that is a prime power.  The modulus of the finite field.
  
Returns:
  - Bool: true or false.
"""
function isgolayperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)
    
	! isprimepower(q) &&  false # we are working in finite fields, so q must be a prime power
    M = q^k
    (isequal(q, 2) && isequal(n, 23) && isequal(d, 7) && isequal(M, 2^12)) && return true
    (isequal(q, 3) && isequal(n, 11) && isequal(d, 5) && isequal(M, 3^6)) && return true
    return false
end

"""
	get_all_words(Σ::Alphabet, q::Integer, ::Val{n}) 		-> Array{Tuple{Symbol}, 1}
	get_all_words(Σ::AbstractArray, q::Integer, n::Integer) -> Array{Tuple{Symbol}, 1}
	get_all_words(Σ::AbstractArray, n::Integer)				-> Array{Tuple{Symbol}, 1}
	get_all_words(q::Integer, n::Integer)					-> Array{Tuple{Symbol}, 1}
	
Get the universe of all codewords of a given alphabet.  The alphabet will be uniquely generated if none is given.
	
Parameters:
  - Σ::AbstractArray: The alphabet allowed.
  - q::Integer: The size of the alphabet.
  - n::Integer: The (fixed) length of the words in the code.
  - d::Integer: The minimum distance between words in the code.
  - 𝒰::AbstractArray: The universe of all codewords of q many letters of block length n.
  
Returns:
  - Array{Tuple{Symbol}, 1}: An array of codewords.  Each codewords is a tuple, and each character in said word is a symbol.
"""
@generated function get_all_words(Σ::Alphabet, q::Integer, ::Val{n}) where n
	quote
		C = Tuple[]
			
		Base.Cartesian.@nloops $n i d -> Σ begin
			wᵢ = Base.Cartesian.@ntuple $n i
			push!(C, wᵢ)
		end
		
		return C
	end
end

get_all_words(Σ::Alphabet, q::Integer, n::Integer) = get_all_words(Σ, q, Val(n))
get_all_words(Σ::Alphabet, n::Integer) = get_all_words(Σ, length(unique(Σ)), n) # if alphabet is given, then q is the length of that alphabet
get_all_words(Σ::AbstractArray, q::Integer, n::Integer) = get_all_words(Alphabet(Σ), q, Val(n))
get_all_words(Σ::AbstractArray, n::Integer) = get_all_words(Alphabet(Σ), length(unique(Σ)), n) # if alphabet is given, then q is the length of that alphabet
get_all_words(q::Integer, n::Integer) = get_all_words(Alphabet(Symbol[gensym() for _ in 1:q]), q, n) # generate symbols if no alphabet is given

"""
	get_codewords_greedy(Σ::AbstractArray, q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray) -> Array{Tuple{Symbol}, 1}
	get_codewords_greedy(Σ::AbstractArray, n::Integer, d::Integer, 𝒰::AbstractArray) -> Array{Tuple{Symbol}, 1}
	get_codewords_greedy(q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray)	-> Array{Tuple{Symbol}, 1}
	get_codewords_greedy(Σ::AbstractArray, q::Integer, n::Integer, d::Integer) -> Array{Tuple{Symbol}, 1}
	get_codewords_greedy(Σ::AbstractArray, n::Integer, d::Integer) -> Array{Tuple{Symbol}, 1}
	get_codewords_greedy(q::Integer, n::Integer, d::Integer) -> Array{Tuple{Symbol}, 1}
	
Search through the universe of all codewords and find a code of block length n and distance d, using the alphabet Σ.  The alphabet will be uniquely generated if none is given.
	
Parameters:
  - Σ::AbstractArray: The alphabet allowed.
  - q::Integer: The size of the alphabet.
  - n::Integer: The (fixed) length of the words in the code.
  - d::Integer: The minimum distance between words in the code.
  - 𝒰::AbstractArray: The universe of all codewords of q many letters of block length n.
  
Returns:
  - Array{Tuple{Symbol}, 1}: An array of codewords.  Each codewords is a tuple, and each character in said word is a symbol.
"""
function get_codewords_greedy(Σ::Alphabet, q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray)
	C = Tuple[]
	
	for wᵢ in 𝒰
		push_if_allowed!(C, wᵢ, d)
	end
	
	return C
end

function get_codewords_greedy(𝒰::UniverseParameters, d::Integer)
	C = Tuple[]
	
	for wᵢ in 𝒰
		push_if_allowed!(C, wᵢ, d)
	end
	
	return C
end

get_codewords_greedy(Σ::Alphabet, q::Integer, n::Integer, d::Integer) =
	get_codewords_greedy(UniverseParameters(Σ, q, n), d)

get_codewords_greedy(Σ::Alphabet, n::Integer, d::Integer) =
	get_codewords_greedy(UniverseParameters(Σ, n), d)

get_codewords_greedy(q::Integer, n::Integer, d::Integer) =
	get_codewords_greedy(UniverseParameters(q, n), d)

get_codewords_greedy(Σ::AbstractArray, q::Integer, n::Integer, d::Integer) =
	get_codewords_greedy(Alphabet(Σ), q, n, d)
	
get_codewords_greedy(Σ::AbstractArray, n::Integer, d::Integer) =
	get_codewords_greedy(UniverseParameters(Σ, n), d)
	
get_codewords_greedy(Σ::AbstractArray, q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray) =
	get_codewords_greedy(Alphabet(Σ), q, n, d, 𝒰)

# if alphabet is given, then q is the length of that alphabet
# get_codewords_greedy(Σ::Alphabet, n::Integer, d::Integer, 𝒰::AbstractArray) =
# 	get_codewords_greedy(Σ, length(unique(Σ)), n, d, 𝒰)
# if the universe of all possible codewords is not given, find it
# get_codewords_greedy(Σ::Alphabet, q::Integer, n::Integer, d::Integer) =
# 	get_codewords_greedy(Σ, q, n, d, get_all_words(Σ, q, n))
# if the universe of all possible codewords is not given, find it and the size of the alphabet
# get_codewords_greedy(Σ::Alphabet, n::Integer, d::Integer) =
# 	get_codewords_greedy(Σ, length(unique(Σ)), n, d, get_all_words(Σ, n))
# if only alphabet size, block length, and distance are given.
# get_codewords_greedy(q::Integer, n::Integer, d::Integer) =
# 	get_codewords_greedy(Alphabet(Symbol[gensym() for _ in 1:q]), q, n, d, get_all_words(q, n))
# if alphabet is given, then q is the length of that alphabet
# get_codewords_greedy(Σ::AbstractArray, n::Integer, d::Integer, 𝒰::AbstractArray) =
# 	get_codewords_greedy(Alphabet(Σ), length(unique(Σ)), n, d, 𝒰)
# generate symbols if no alphabet is given
# get_codewords_greedy(q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray) =
# 	get_codewords_greedy(Alphabet(Symbol[gensym() for _ in 1:q]), q, n, d, 𝒰)
# if the universe of all possible codewords is not given, find it
# get_codewords_greedy(Σ::AbstractArray, q::Integer, n::Integer, d::Integer) =
# 	get_codewords_greedy(Alphabet(Σ), q, n, d, get_all_words(Σ, q, n))
# if the universe of all possible codewords is not given, find it and the size of the alphabet
# get_codewords_greedy(Σ::AbstractArray, n::Integer, d::Integer) =
# 	get_codewords_greedy(Alphabet(Σ), length(unique(Σ)), n, d, get_all_words(Σ, n))
# if only alphabet size, block length, and distance are given.
# get_codewords_greedy(q::Integer, n::Integer, d::Integer) =
# 	get_codewords_greedy(Symbol[gensym() for _ in 1:q], q, n, d, get_all_words(q, n))

"""
	get_codewords_random(Σ::AbstractArray, q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray) -> Array{Tuple{Symbol}, 1}
	get_codewords_random(Σ::AbstractArray, n::Integer, d::Integer, 𝒰::AbstractArray) -> Array{Tuple{Symbol}, 1}
	get_codewords_random(q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray)	-> Array{Tuple{Symbol}, 1}
	get_codewords_random(Σ::AbstractArray, q::Integer, n::Integer, d::Integer) -> Array{Tuple{Symbol}, 1}
	get_codewords_random(Σ::AbstractArray, n::Integer, d::Integer) -> Array{Tuple{Symbol}, 1}
	get_codewords_random(q::Integer, n::Integer, d::Integer) -> Array{Tuple{Symbol}, 1}
	
Search through the universe of all codewords at random and find a code of block length n and distance d, using the alphabet Σ.  The alphabet will be uniquely generated if none is given.
	
Parameters:
  - Σ::AbstractArray: The alphabet allowed.
  - q::Integer: The size of the alphabet.
  - n::Integer: The (fixed) length of the words in the code.
  - d::Integer: The minimum distance between words in the code.
  - 𝒰::AbstractArray: The universe of all codewords of q many letters of block length n.
  
Returns:
  - Array{Tuple{Symbol}, 1}: An array of codewords.  Each codewords is a tuple, and each character in said word is a symbol.
"""
function get_codewords_random(Σ::AbstractArray, q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray)
	C = Tuple[]
	Σ = ensure_symbolic(Σ)
	Σ = unique(Σ)
	𝒰′ = copy(𝒰)
	
	while ! isempty(𝒰′)
		wᵢ = rand(𝒰′)
		
		push_if_allowed!(C, wᵢ, d)
		deleteat!(𝒰′, findfirst(x -> isequal(x, wᵢ), 𝒰′))
	end
	
	return C
end

# if alphabet is given, then q is the length of that alphabet
get_codewords_random(Σ::AbstractArray, n::Integer, d::Integer, 𝒰::AbstractArray) =
	get_codewords_random(Σ, length(unique(Σ)), n, d, 𝒰)
# generate symbols if no alphabet is given
get_codewords_random(q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray) =
	get_codewords_random(Symbol[gensym() for _ in 1:q], q, n, d, 𝒰)
# if the universe of all possible codewords is not given, find it
get_codewords_random(Σ::AbstractArray, q::Integer, n::Integer, d::Integer) =
	get_codewords_random(Σ, q, n, d, get_all_words(Σ, q, n))
# if the universe of all possible codewords is not given, find it and the size of the alphabet
get_codewords_random(Σ::AbstractArray, n::Integer, d::Integer) =
	get_codewords_random(Σ, length(unique(Σ)), n, d, get_all_words(Σ, n))
# if only alphabet size, block length, and distance are given.
get_codewords_random(q::Integer, n::Integer, d::Integer) =
	get_codewords_random(Symbol[gensym() for _ in 1:q], q, n, d, get_all_words(q, n))

"""
	get_codewords(Σ::AbstractArray, q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray; m::Integer=10) -> Array{Tuple{Symbol}, 1}
	get_codewords(Σ::AbstractArray, n::Integer, d::Integer, 𝒰::AbstractArray; m::Integer=10) -> Array{Tuple{Symbol}, 1}
	get_codewords(q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray; m::Integer=10) -> Array{Tuple{Symbol}, 1}
	get_codewords(Σ::AbstractArray, q::Integer, n::Integer, d::Integer; m::Integer=10) -> Array{Tuple{Symbol}, 1}
	get_codewords(Σ::AbstractArray, n::Integer, d::Integer; m::Integer=10) -> Array{Tuple{Symbol}, 1}
	get_codewords(q::Integer, n::Integer, d::Integer; m::Integer=10) -> Array{Tuple{Symbol}, 1}

Use function `get_codewords_random` m many times, and `get_codewords_greedy`.  Return the code with the greatest number of words.  The alphabet will be uniquely generated if none is given.  You can omit Σ and 𝒰.  You can omit q if Σ is given.
	
Parameters:
  - Σ::AbstractArray: The alphabet allowed.
  - q::Integer: The size of the alphabet.
  - n::Integer: The (fixed) length of the words in the code.
  - d::Integer: The minimum distance between words in the code.
  - 𝒰::AbstractArray: The universe of all codewords of q many letters of block length n.
  - m::Integer (kwarg): Try a random code m many times.
  
Returns:
  - Array{Tuple{Symbol}, 1}: An array of codewords.  Each codewords is a tuple, and each character in said word is a symbol.
"""
function get_codewords(Σ::AbstractArray, q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray; m::Integer=10)
	code_size = 0
	C = Tuple[]
	Σ = ensure_symbolic(Σ)
	Σ = unique(Σ)
	

	for _ in 1:m
		random_code = get_codewords_random(Σ, q, n, d, 𝒰)
		random_size = length(random_code)
		if random_size > code_size
			code_size = random_size
			C = random_code
		end
	end
	
	greedy_code = get_codewords_greedy(Σ, q, n, d, 𝒰)
	greedy_size = length(greedy_code)
	if greedy_size > code_size
		code_size = greedy_size
		C = greedy_code
	end
	
	return C
end

# if alphabet is given, then q is the length of that alphabet
get_codewords(Σ::AbstractArray, n::Integer, d::Integer, 𝒰::AbstractArray; m::Integer=10) =
	get_codewords(Σ, length(unique(Σ)), n, d, 𝒰, m = m)
# generate symbols if no alphabet is given
get_codewords(q::Integer, n::Integer, d::Integer, 𝒰::AbstractArray; m::Integer=10) =
	get_codewords(Symbol[gensym() for _ in 1:q], q, n, d, 𝒰, m = m)
# if the universe of all possible codewords is not given, find it
get_codewords(Σ::AbstractArray, q::Integer, n::Integer, d::Integer; m::Integer=10) =
	get_codewords(Σ, q, n, d, get_all_words(Σ, q, n), m = m)
# if the universe of all possible codewords is not given, find it and the size of the alphabet
get_codewords(Σ::AbstractArray, n::Integer, d::Integer; m::Integer=10) =
	get_codewords(Σ, length(unique(Σ)), n, d, get_all_words(Σ, n), m = m)
# if only alphabet size, block length, and distance are given.
get_codewords(q::Integer, n::Integer, d::Integer; m::Integer=10) =
	get_codewords(Symbol[gensym() for _ in 1:q], q, n, d, get_all_words(q, n), m = m)

"""
	get_codewords(G::AbstractArray, m::Integer) -> Array{Tuple{Symbol}, 1}

Get codewords of a code from the generating matrix under a finite field of modulo m.  Precisely, computes all linear combinations of the rows of the generating matrix.
	
Parameters:
  - G::AbstractArray: A matrix of integers which generates the code.
  - m::Integer: The bounds of the finite field (i.e., the molulus you wish to work in).
  
Returns:
  - Array{Tuple{Symbol}, 1}: An array of codewords.  Each codewords is a tuple, and each character in said word is a symbol.
"""
function get_codewords(G::AbstractArray, m::Integer)
	codewords = Vector()
	rows = Vector(undef, size(G, 2))
	
	for i in 1:size(G, 1)
		rows[i] = [G[i, j] for j in 1:size(G, 2)]
	end
	
	for c in Base.Iterators.product([0:m-1 for i in 1:size(G, 1)]...)
		word = Ref(c[1]) .* rows[1]
		
		for i in 2:size(G, 1)
			word = mod.(word .+ (Ref(c[i]) .* rows[i]), m)
		end
		
		push!(codewords, word)
	end
	
	return codewords
end

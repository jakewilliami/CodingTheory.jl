#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

include("distance.jl")
include("primes.jl")

"""
```julia
push_if_allowed!(C::AbstractArray{T}, w::T, d::Int)
```

Takes in an array and a word.  As long as the word does not mean that the distance is smaller than d, we add w to the array.  If we are successful in doing this, return true.  Otherwise, return false.  *This is a mutating function.*
"""
function push_if_allowed!(C::AbstractArray{T}, w::T, d::Int) where T <: AbstractWord
	isempty(C) && (push!(C, w); return true)
	
	for c in C
		if hamming_distance(c, w) < d
			return false
		end
	end
	
	push!(C, w)
	return true
end

"""
```julia
push_if_allowed!(C::AbstractArray{T}, C′::AbstractArray{T}, w::T, d::Int)
```

Takes in two arrays, A and B.  If w is allowed in C given distance d, push to C′.  If we are successful in doing this, return true.  Otherwise, return false.  *This is a mutating function.*
"""
function push_if_allowed!(C::AbstractVector{T}, C′::AbstractVector{T}, w::T, d::Int) where T <: AbstractWord
	isempty(C) && (push!(C′, w); return true)
	
	for c in C
		if hamming_distance(c, w) < d
			return false
		end
	end
	
	push!(C′, w)
	return true
end

"""
```julia
replace_if_allowed!(C::AbstractArray, d::Int, w, w′) -> Bool
```

Takes in an array and a word.  As long as the word does not mean that the distance is smaller than d, we replace a with b in the array.  Replaces and returns true if allowed; otherwise returns false.  *This is a mutating function.*
"""
function replace_if_allowed!(C::AbstractVector{T}, d::Int, ws::Pair) where T <: AbstractWord
	w, w′ = ws
	
	for c in C
		if hamming_distance(c, w′) < d
			return false
		end
	end
	
	replace!(C, ws)
	return true
end
function replace_if_allowed!(C::AbstractVector{T}, d::Int, w::T, w′::T) where T <: AbstractWord
	return replace_if_allowed!(C, d, Pair(w, w′))
end

_mutate_codeword(w::Word{N, T}, i::Int, a::T) where {T, N} =
	setindex!(w, a, i)
	
"""
```julia
mutate_codeword(w::NonStaticAbstractWord{N, T}, n::Int, i::Int, a::T) where {T, N} -> MVector{N, T}
mutate_codeword(w::Word{N, T}, n::Int, i::Int, a::T) where {T, N} -> MVector{N, T}
```

Mutates the word w, which is an `MVector` of length N, changing its iᵗʰ index to a.
"""
mutate_codeword(w::NonStaticAbstractWord{N, T}, i::Int, a::T) where {T, N} =
	_mutate_codeword(Word{N, T}(w), i, a)
mutate_codeword(w::Word{N, T}, n::Int, i::Int, a::T) where {T, N} =
	_mutate_codeword(w, i, a)

"""
```julia
get_all_words(Σ::Alphabet{N}, q::Int, n::Int) -> Codewords{M}
get_all_words(Σ::Alphabet{N}, n::Int) -> Codewords{M}
get_all_words(Σ::AbstractArray, q::Int, n::Int) -> Codewords{M}
get_all_words(Σ::AbstractArray, n::Int) -> Codewords{M}
get_all_words(q::Int, n::Int) -> Codewords{M}
```
	
Get the universe of *all* codewords of a given alphabet.  The alphabet will be uniquely generated if none is given.
	
Parameters:
  - `Σ::AbstractArray`: The alphabet allowed.
  - `q::Int`: The size of the alphabet.
  - `n::Int`: The (fixed) length of the words in the code.
  - `d::Int`: The minimum distance between words in the code.
  - `𝒰::AbstractArray`: The universe of all codewords of q many letters of block length n.
  
Returns:
  - `Codewords{M}`: An array of codewords, each of length `M`.  Each codewords is a tuple, and each character in said word is a symbol.

---

### Examples

```julia
julia> get_all_words(2, 2) # all words of block length 2 using 2 unique symbols
2×2 Array{Tuple{Symbol,Symbol},2}:
 (Symbol("##254"), Symbol("##254"))  (Symbol("##254"), Symbol("##253"))
 (Symbol("##253"), Symbol("##254"))  (Symbol("##253"), Symbol("##253"))
```
"""
get_all_words(Σ::Alphabet{N}, q::Int, n::Int) where {N} =
	collect(CodeUniverseIterator(UniverseParameters(Σ, q, n)))
get_all_words(Σ::Alphabet{N}, n::Int) where {N} =
	get_all_words(Σ, length(Set(Σ)), n) # if alphabet is given, then q is the length of that alphabet
get_all_words(Σ::AbstractArray, n::Int) =
	get_all_words(Alphabet(Σ), length(Set(Σ)), n) # if alphabet is given, then q is the length of that alphabet
get_all_words(q::Int, n::Int) =
	get_all_words(genalphabet(q), q, n) # generate symbols if no alphabet is given

"""
```julia
get_codewords_greedy(𝒰::UniverseParameters, d::Int) -> Codewords{M}
get_codewords_greedy(Σ::Alphabet{N}, q::Int, n::Int, d::Int) -> Codewords{M}
get_codewords_greedy(Σ::Alphabet{N}, n::Int, d::Int) -> Codewords{M}
get_codewords_greedy(q::Int, n::Int, d::Int) -> Codewords{M}
get_codewords_greedy(Σ::AbstractArray, q::Int, n::Int, d::Int) -> Codewords{M}
get_codewords_greedy(Σ::AbstractArray, n::Int, d::Int) -> Codewords{M}
get_codewords_greedy(Σ::AbstractArray, q::Int, n::Int, d::Int, 𝒰::AbstractArray) -> Codewords{M}
```
	
Search through the universe of all codewords and find a code of block length n and distance d, using the alphabet Σ.  The alphabet will be uniquely generated if none is given.  This uses a greedy algorithm, simply iterating through all words (see above) and choosing them if they fit in the code.  In some cases the greedy algorithm is the best, but in others it is very much not.
	
Parameters:
  - `𝒰::UniverseParameters`: The parameters of the universe of all codewords of q many letters of block length n.
  - `Σ::AbstractArray`: The alphabet allowed.
  - `q::Int`: The size of the alphabet.
  - `n::Int`: The (fixed) length of the words in the code.
  - `d::Int`: The minimum distance between words in the code.
  
Returns:
  - `Codewords{M}`: An array of codewords, each of length `M`.  Each codewords is a tuple, and each character in said word is a symbol.
"""
function get_codewords_greedy(𝒰::UniverseParameters, d::Int)
	C = eltype(𝒰)[]
	
	for wᵢ in CodeUniverseIterator(𝒰)
		push_if_allowed!(C, wᵢ, d)
	end
	
	return C
end
get_codewords_greedy(Σ::Alphabet{N}, q::Int, n::Int, d::Int) where {N} =
	get_codewords_greedy(UniverseParameters(Σ, q, n), d)
get_codewords_greedy(Σ::Alphabet{N}, n::Int, d::Int) where {N} =
	get_codewords_greedy(UniverseParameters(Σ, n), d)
get_codewords_greedy(q::Int, n::Int, d::Int) =
	get_codewords_greedy(UniverseParameters(q, n), d)
get_codewords_greedy(Σ::AbstractArray, q::Int, n::Int, d::Int) =
	get_codewords_greedy(Alphabet(Σ), q, n, d)
get_codewords_greedy(Σ::AbstractArray, n::Int, d::Int) =
	get_codewords_greedy(UniverseParameters(Σ, n), d)
get_codewords_greedy(Σ::AbstractArray, q::Int, n::Int, d::Int, 𝒰::AbstractArray) =
	get_codewords_greedy(Alphabet(Σ), q, n, d, 𝒰)
	

argmaxminima(A::AbstractArray; dims::Int) = getindex(argmin(A, dims = dims), argmax(argmin(A, dims = dims)))
maxminima(A::AbstractArray; dims::Int) = getindex(minimum(A, dims = dims), maximum(minimum(A, dims = dims)))
argminmaxima(A::AbstractArray; dims::Int) = getindex(argmax(A, dims = dims), argmin(argmax(A, dims = dims)))
minmaxima(A::AbstractArray; dims::Int) = getindex(maximum(A, dims = dims), minimum(maximum(A, dims = dims)))

"""
```julia
get_codewords_random(𝒰::UniverseParameters, d::Int; m::Int = 1000) -> Codewords{M}
get_codewords_random(Σ::Alphabet{N}, q::Int, n::Int, d::Int; m::Int=1000) -> Codewords{M}
get_codewords_random(Σ::Alphabet{N}, n::Int, d::Int; m::Int=1000)	-> Codewords{M}
get_codewords_random(q::Int, n::Int, d::Int; m::Int=1000) -> Codewords{M}
get_codewords_random(Σ::AbstractArray, q::Int, n::Int, d::Int; m::Int=1000) -> Codewords{M}
get_codewords_random(Σ::AbstractArray, n::Int, d::Int; m::Int=1000) -> Codewords{M}
get_codewords_random(Σ::AbstractArray, q::Int, n::Int, d::Int, 𝒰::AbstractArray; m::Int=1000) -> Codewords{M}
```

Search through the universe of all codewords at random and find a code of block length n and distance d, using the alphabet Σ.  The alphabet will be uniquely generated if none is given.  This is a cleverer algorithm than the greedy algorithm.  Increasing the `m` keyword argument arbitrarily _should_ produce a maximal code, as for each codeword it chooses, it collects a list of `m` many random words, and chooses the best one from that intermediate list.

Parameters:
  - `Σ::AbstractArray`: The alphabet allowed.
  - `q::Int`: The size of the alphabet.
  - `n::Int`: The (fixed) length of the words in the code.
  - `d::Int`: The minimum distance between words in the code.
  - `𝒰::AbstractArray`: The universe of all codewords of q many letters of block length n.

Returns:
  - `Codewords{M}`: An array of codewords, each of length `M`.  Each codewords is a tuple, and each character in said word is a symbol.
"""
function get_codewords_random(𝒰::UniverseParameters, d::Int; m::Int = 1000)
	C = eltype(𝒰)[]
	
	starting_word = rand(𝒰) # get a random word in the code start
	push!(C, starting_word)
	
	for _ in 1:length(𝒰)
		C′ = eltype(𝒰)[]
		# while length(C′) < m
		for _ in 1:m
			wᵣ = rand(𝒰)
			push_if_allowed!(C, C′, wᵣ, d) # if allowed in C, push to C′
			# wᵣ ∉ C′ && push_if_allowed!(C, C′, wᵣ, d)
		end
		isempty(C′) && break
		# [push_if_allowed!(C, C′, w, d) for _ in 1:m]
		distances = Int[hamming_distance(wᵢ, wⱼ) for wᵢ in C, wⱼ in C′]
		best_word = getindex(C′, getindex(argmaxminima(distances, dims = 1), 2))
		push!(C, best_word)
	end
	
	return C
end
get_codewords_random(Σ::Alphabet{N}, q::Int, n::Int, d::Int; m::Int=1000) where {N} =
	get_codewords_random(UniverseParameters(Σ, q, n), d, m=m)
get_codewords_random(Σ::Alphabet{N}, n::Int, d::Int; m::Int=1000) where {N} =
	get_codewords_random(UniverseParameters(Σ, n), d, m=m)
get_codewords_random(q::Int, n::Int, d::Int; m::Int=1000) =
	get_codewords_random(UniverseParameters(q, n), d, m=m)
get_codewords_random(Σ::AbstractArray, q::Int, n::Int, d::Int; m::Int=1000) =
	get_codewords_random(Alphabet(Σ), q, n, d, m=m)
get_codewords_random(Σ::AbstractArray, n::Int, d::Int; m::Int=1000) =
	get_codewords_random(UniverseParameters(Σ, n), d, m=m)
get_codewords_random(Σ::AbstractArray, q::Int, n::Int, d::Int, 𝒰::AbstractArray; m::Int=1000) =
	get_codewords_random(Alphabet(Σ), q, n, d, 𝒰, m=m)


# using Mmap
# function get_codewords_random_mmap(mmappath::AbstractString, 𝒰::UniverseParameters, d::Int)
# 	io = open(mmappath, "r+") # allow read and write
# 	# write(s, size(A,2))
# 	# read(s, Int)
# 	# close(s)
# 	# B = Mmap.mmap(io, BitArray, (25,30000))
# 	# Mmap.sync!(B);
# 	# close(io);
# 	# rm("mmap.bin")
# 	starting_word = rand(𝒰) # get a random word in the code start
# 	push!(C, starting_word)
#
# 	for _ in 1:length(𝒰)
# 		C′ = Tuple[]
# 		for _ in 1:m
# 			push_if_allowed!(C, C′, rand(𝒰), d) # if allowed in C, push to C′
# 		end
# 		isempty(C′) && break
# 		# [push_if_allowed!(C, C′, w, d) for _ in 1:m]
# 		distances = [hamming_distance(wᵢ, wⱼ) for wᵢ in C, wⱼ in C′]
# 		best_word = getindex(C′, getindex(argmaxminima(distances, dims = 1), 2))
# 		push!(C, best_word)
# 	end
#
# 	return C
# end

# get_codewords_random(𝒰::UniverseParameters, d::Int) = get_codewords_random(joinpath(tempdir(), "mmap.bin"), 𝒰, d)

"""
```julia
get_codewords(𝒰::UniverseParameters, d::Int; m::Int=10, m_random::Int = 1000) -> Codewords{M}
get_codewords(Σ::Alphabet{N}, q::Int, n::Int, d::Int; m::Int=10, m_random::Int=1000) -> Codewords{M}
get_codewords(q::Int, n::Int, d::Int; m::Int=10, m_random::Int=1000) -> Codewords{M}
get_codewords(q::Int, n::Int, d::Int; m::Int=10, m_random::Int=1000) -> Codewords{M}
get_codewords(Σ::AbstractArray, q::Int, n::Int, d::Int; m::Int=10, m_random::Int=1000) -> Codewords{M}
get_codewords(Σ::AbstractArray, n::Int, d::Int; m::Int=10, m_random::Int=1000) -> Codewords{M}
get_codewords(Σ::AbstractArray, q::Int, n::Int, d::Int, 𝒰::AbstractArray; m::Int=10, m_random::Int=1000) -> Codewords{M}
```

Use function `get_codewords_random` `m` many times (with `get_codewords_random(..., m = m_random)`), and `get_codewords_greedy`.  Return the code with the greatest number of words.  The alphabet will be uniquely generated if none is given.  You can omit Σ and 𝒰.  You can omit q if Σ is given.
	
Parameters:
  - `Σ::AbstractArray`: The alphabet allowed.
  - `q::Int`: The size of the alphabet.
  - `n::Int`: The (fixed) length of the words in the code.
  - `d::Int`: The minimum distance between words in the code.
  - `𝒰::AbstractArray`: The universe of all codewords of q many letters of block length n.
  - `m::Int` (kwarg): Try a random code m many times.
  - `m_random::Int` (kwarg): The number of possible words `get_codewords_random` chooses from for _each_ word it selects.
  
Returns:
  - `Codewords{M}`: An array of codewords, each of length `M`.  Each codewords is a tuple, and each character in said word is a symbol.

!!! note

    If you are looking for a _maximal_ code, this is likely the function you need.  Increasing `m` and `m_random` arbitrarily should ensure a maximal code—_however_, that computing power/time in not always possible, as it requires a lot of RAM to store certain codes in memeory.  Efforts are being made to make this process better by using memory-mapped filed instead of storing codewords in RAM, but this will make it much slower as well.  Help with this would be much appreciated.

---

### Examples

```julia
julia> get_codewords(["a", "b", "c"], 3, 2) # get codewords of block length 3 with distance 2.  Once again, are symbols for uniqueness
9-element Array{Tuple{Symbol,Symbol,Symbol},1}:
 (:a, :b, :a)
 (:c, :a, :b)
 (:b, :c, :c)
 (:b, :a, :a)
 (:c, :b, :c)
 (:a, :a, :c)
 (:a, :c, :b)
 (:c, :c, :a)
 (:b, :b, :b)
```
"""
function get_codewords(𝒰::UniverseParameters, d::Int; m::Int=10, m_random::Int = 1000)
	code_size = 0
	C = eltype(𝒰)[]
	

	for _ in 1:m
		random_code = get_codewords_random(𝒰, d; m = m_random)
		random_size = length(random_code)
		if random_size > code_size
			code_size = random_size
			C = random_code
		end
	end
	
	greedy_code = get_codewords_greedy(𝒰, d)
	greedy_size = length(greedy_code)
	if greedy_size > code_size
		code_size = greedy_size
		C = greedy_code
	end
	
	return C
end
get_codewords(Σ::Alphabet{N}, q::Int, n::Int, d::Int; m::Int=10, m_random::Int=1000) where {N} =
	get_codewords(UniverseParameters(Σ, q, n), d, m=m, m_random=m_random)
get_codewords(Σ::Alphabet{N}, n::Int, d::Int; m::Int=10, m_random::Int=1000) where {N} =
	get_codewords(UniverseParameters(Σ, n), d, m=m, m_random=m_random)
get_codewords(q::Int, n::Int, d::Int; m::Int=10, m_random::Int=1000) =
	get_codewords(UniverseParameters(q, n), d, m=m, m_random=m_random)
get_codewords(Σ::AbstractArray, q::Int, n::Int, d::Int; m::Int=10, m_random::Int=1000) =
	get_codewords(Alphabet(Σ), q, n, d, m=m, m_random=m_random)
get_codewords(Σ::AbstractArray, n::Int, d::Int; m::Int=10, m_random::Int=1000) =
	get_codewords(UniverseParameters(Σ, n), d, m=m, m_random=m_random)
get_codewords(Σ::AbstractArray, q::Int, n::Int, d::Int, 𝒰::AbstractArray; m::Int=10, m_random::Int=1000) =
	get_codewords(Alphabet(Σ), q, n, d, 𝒰, m=m, m_random=m_random)

"""
```julia
get_codewords(G::AbstractArray, m::Int) -> Codewords{M}
```

Get codewords of a code from the _generating matrix_ under a finite field of modulo `m`.  Precisely, computes all linear combinations of the rows of the generating matrix.
	
Parameters:
  - `G::AbstractArray`: A matrix of Ints which generates the code.
  - `m::Int`: The bounds of the finite field (i.e., the molulus you wish to work in).
  
Returns:
  - `Codewords{M}`: An array of codewords, each of length `M`.  Each codewords is a tuple, and each character in said word is a symbol.
"""
function get_codewords(G::AbstractArray, m::Int)
	codewords = Vector()
	rows = Vector(undef, size(G, 2))
	
	for i in 1:size(G, 1)
		rows[i] = [G[i, j] for j in 1:size(G, 2)]
	end
	
	for c in Base.Iterators.product([0:m-1 for _ in 1:size(G, 1)]...)
		word = Ref(c[1]) .* rows[1]
		
		for i in 2:size(G, 1)
			word = mod.(word .+ (Ref(c[i]) .* rows[i]), m)
		end
		
		push!(codewords, word)
	end
	
	return codewords
end

# function obtain_maximal_code(𝒰::UniverseParameters, d::Int)
# 	adj_matrix = Matrix{Int8}(undef, length(𝒰), length(𝒰))
#
# 	for u in CodeUniverseIterator(𝒰), u′ in CodeUniverseIterator(𝒰)
# 		distance = hamming_distance(u, u′)
# 	end
# end

#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "distance.jl"))
include(joinpath(dirname(@__FILE__), "utils.jl"))
include(joinpath(dirname(@__FILE__), "primes.jl"))
include(joinpath(dirname(@__FILE__), "rref.jl"))

using LinearAlgebra: I

struct Alphabet
    Σ::Union{AbstractArray, AbstractString}
    
    function Alphabet(Σ::AbstractArray)
        try
            Σ = parse.(Int, Σ)
        catch
        end
        
        new(Σ)
    end # end constructor function
    
    function Alphabet(Σ::AbstractString)
        Σ = collect(Σ)
        
        try
            Σ = parse.(Int, Σ)
        catch
            Σ = string.(Σ)
        end
        
        new(Σ)
    end # end constructor function
end # end struct

struct Messages
    ℳ::AbstractArray
    length::Integer
    
    function Messages(ℳ::AbstractArray)
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(ℳ) || throw(error("$(message_length_error)"))
        
        length = length(ℳ[rand(1:length(ℳ))]) # choose arbitrary message in the list of messages
        
        new(ℳ, length)
    end # end constructor function
end

#=
n is the word length
q is the number of symbols in the code
M is the size/number of elements in the code
=#
rate(q::Integer, M::Integer, n::Integer) = log(q, M) / n

struct Rounding end
no_round = Rounding()

__spheres(q::Integer, n::Integer, r::Integer) = sum([((big(q) - 1)^i) * binomial(big(n), big(i)) for i in 0:r])
__sphere_bound(round_func::Function, q::Integer, n::Integer, d::Integer) = round_func((big(q)^n) / __spheres(q, n, d))
sphere_covering_bound(q::Integer, n::Integer, d::Integer)::Integer = __sphere_bound(ceil, q, n, d - 1)
sphere_packing_bound(q::Integer, n::Integer, d::Integer)::Integer = __sphere_bound(floor, q, n, Int(floor((d - 1) / 2)))
sphere_packing_bound(q::Integer, n::Integer, d::Integer, ::Rounding)::Real = __sphere_bound(identity, q, n, Int(floor((d - 1) / 2)))
hamming_bound(q::Integer, n::Integer, d::Integer)::Integer = sphere_packing_bound(q, n, d)
hamming_bound(q::Integer, n::Integer, d::Integer, ::Rounding)::Real = sphere_packing_bound(q, n, d, no_round)
singleton_bound(q::Number, n::Number, d::Number)::Real = float(big(q))^(big(n) - big(d) + 1)

function construct_ham_matrix(r::Integer, q::Integer)::Matrix
    ncols = Int(floor((q^r - 1) / (q - 1)))
    M = Matrix{Integer}(undef, r, ncols)
    
    for i in 1:ncols
        M[:,i] = reverse(digits(parse(Int, string(i, base = q)), pad = r), dims = 1)
    end
    
    return M
end

#=
n is the word length
k is the dimension of the code (i.e., )
d is the distance of the code
If C ⊆ 𝔽_q^n, then C is over the field modulo q
=#
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

function ishammingperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)::Bool
    
    isprimepower(q) || return false
    
    M = q^k
    r = log(ℯ, ((n * log(ℯ, 1)) / (log(ℯ, 2))) + 1) / log(ℯ, 2)
        
    if isequal(n, (q^(r - 1)) / (q - 1)) && isequal(d, 3) && isequal(M, q^(((q^r - 1) / (q - 1)) - r))
        return true
    end
    
    return false
end

function ishammingbound(q::Integer, n::Integer, d::Integer)::Bool
	! isprimepower(q) && return false
	d ≠ 3 && return false
	
	r = 1
	while ((q^r - 1) / (q - 1)) < n
		r = r + 1
	end
	
	return isequal(((q^r - 1) / (q - 1)), n) ? true : false
end

function isgolayperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)::Bool
    
    M = q^k
    
    if isequal(q, 2) && isequal(n, 23) && isequal(d, 7) && isequal(M, 2^12)
        return true
    end
    
    if isequal(q, 3) && isequal(n, 11) && isequal(d, 5) && isequal(M, 3^6)
        return true
    end
    
    return false
end

#=
Finds all possible combinations of words of length n using q symbols from alphabet Σ
=#
@generated function get_all_words(Σ::AbstractArray, q::Integer, ::Val{n}) where n
	quote
		C = Tuple[]
		Σ = __ensure_symbolic(Σ)
		Σ = unique(Σ)
			
		Base.Cartesian.@nloops $n i d -> Σ begin
			wᵢ = Base.Cartesian.@ntuple $n i
			push!(C, wᵢ)
		end
		
		return C
	end
end

get_all_words(Σ::AbstractArray, q::Integer, n::Integer) = get_all_words(Σ, q, Val(n))
get_all_words(Σ::AbstractArray, n::Integer) = get_all_words(Σ, length(unique(Σ)), n) # if alphabet is given, then q is the length of that alphabet
get_all_words(q::Integer, n::Integer) = get_all_words(Symbol[gensym() for _ in 1:q], q, n) # generate symbols if no alphabet is given

#=
Get a list of messages with q letters in the alphabet, and word size of n
=#
function get_codewords_greedy(Σ::AbstractArray, q::Integer, n::Integer, d::Integer)
	C = Tuple[]
	Σ = __ensure_symbolic(Σ)
	Σ = unique(Σ)
	all_codewords = get_all_words(Σ, q, n)
	
	for wᵢ in all_codewords
		__push_if_allowed!(C, wᵢ, d)
	end
	
	return C
end

get_codewords_greedy(Σ::AbstractArray, n::Integer, d::Integer) = get_codewords_greedy(Σ, length(unique(Σ)), n, d) # if alphabet is given, then q is the length of that alphabet
get_codewords_greedy(q::Integer, n::Integer, d::Integer) = get_codewords_greedy(Symbol[gensym() for _ in 1:q], q, n, d) # generate symbols if no alphabet is given

function get_codewords_random(Σ::AbstractArray, q::Integer, n::Integer, d::Integer)
	C = Tuple[]
	Σ = __ensure_symbolic(Σ)
	Σ = unique(Σ)
	all_codewords = get_all_words(Σ, n)
	
	while ! isempty(all_codewords)
		wᵢ = rand(all_codewords)
		
		__push_if_allowed!(C, wᵢ, d)
		deleteat!(all_codewords, findall(x -> isequal(x, wᵢ), all_codewords))
	end
	
	return C
end

get_codewords_random(Σ::AbstractArray, n::Integer, d::Integer) = get_codewords_random(Σ, length(unique(Σ)), n, d) # if alphabet is given, then q is the length of that alphabet
get_codewords_random(q::Integer, n::Integer, d::Integer) = get_codewords_random(Symbol[gensym() for _ in 1:q], q, n, d) # generate symbols if no alphabet is given

function get_codewords(Σ::AbstractArray, q::Integer, n::Integer, d::Integer; m::Integer=10)
	size = 0
	C = Tuple[]
	Σ = __ensure_symbolic(Σ)
	Σ = unique(Σ)
		
	for _ in 1:m
		random_code = get_codewords_random(Σ, q, n, d)
		random_size = length(random_code)
		if random_size > size
			size = random_size
			C = random_code
		end
	end
	
	greedy_code = get_codewords_greedy(Σ, q, n, d)
	greedy_size = length(greedy_code)
	if greedy_size > size
		size = greedy_size
		C = greedy_code
	end
	
	return C
end

get_codewords(Σ::AbstractArray, n::Integer,  d::Integer; m::Integer=10) = get_codewords(Σ, length(unique(Σ)), n, d, m=m) # if alphabet is given, then q is the length of that alphabet
get_codewords(q::Integer, n::Integer,  d::Integer; m::Integer=10) = get_codewords(Symbol[gensym() for _ in 1:q], q, n, d, m=m) # generate symbols if no alphabet is given

#=
Get codewords by computing all linear combinations of the rows of the matrix, under modulo
=#
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

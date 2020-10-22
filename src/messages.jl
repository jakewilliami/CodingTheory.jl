#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "distance.jl"))
include(joinpath(dirname(@__FILE__), "utils.jl"))
include(joinpath(dirname(@__FILE__), "primes.jl"))

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
    M = Matrix{Int}(undef, r, ncols)
    
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
Get a list of messages with q letters in the alphabet, and word size of n

E.g., if (q, n, d) = (3, 2, 1), get_alphabet returns:
9-element Array{String,1}:
 "11"
 "21"
 "31"
 "12"
 "22"
 "32"
 "13"
 "23"
 "33"
And expands to

=#
@generated function get_codewords(Σ::AbstractArray, q::Integer, ::Val{n}, d::Integer)::Array{NTuple{Symbol, N}, 1} where n
	quote
		C = Tuple[]
		
		if eltype(Σ) isa Symbol
		else
			Σ = __deepsym(Σ)
		end
			
		Base.Cartesian.@nloops $n i d -> Σ begin
			d > n && continue
			wᵢ = Base.Cartesian.@ntuple $n i
			
			if isempty(C) || code_distance(C, wᵢ) ≥ d
				push!(C, wᵢ)
			end
		end
		
		return C
	end
end

get_codewords(Σ::AbstractArray, q::Integer, n::Integer, d::Integer) = get_codewords(Σ, q, Val(n), d::Integer)
get_codewords(q::Integer, n::Integer, d::Integer) = get_codewords(Symbol[gensym() for _ in 1:q], q, n, d)

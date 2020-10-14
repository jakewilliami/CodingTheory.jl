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

function perfect_search()
    q = 6 # a non prime power
    count = 0
    n, k, d = 0, 0, 0

    for n in 3:100, k in 3:100, d in 3:100
        # println(n,k,d,q)
        if isperfect(n, k, d, q) #&& ! ishammingperfect(n, k, d, q) && ! isgolayperfect(n, k, d, q)
            # return n, k, d, q
            println("[$(n), $(k), $(d)]$(q)-code")
        end
    end

    return nothing
end


using CSV, DataFrames, StatsPlots
# You might want to use your code to search for parameters where q^n divided by sphere_packing_bound(q,n,d) is an integer. I think this is essentially how Golay started his search - it would be interesting to replicate some of that search.
function integer_search()::Array{Array{Int, 1}}
    A = []
    upper_bound = 10
    increment_bound = 10
    stop_at = 1_000_000
    
    while true
        for q in 1:upper_bound, n in 1:upper_bound, d in 1:upper_bound
            __arelessthan(upper_bound - increment_bound, q, n, d) && continue
            # 1 ∈ (q, n, d) && continue
            # 2 ∈ (q, n, d) && continue
            if isprimepower(q)
                if __aredistinct(q, n, d)
                # if ! __allequal(q, n, d)
                    hb = hamming_bound(q, n, d, no_round)
                    if isinteger(hb)
                        if ! isone(hb)
                            push!(A, [q, n, d, hb])
                            isequal(length(A), stop_at) && return A
                        end
                    end
                end
            end
        end
        upper_bound += increment_bound
    end
end

function make_csv()
    A = integer_search()
    D = DataFrame(q = Number[], n = Number[], d = Number[], hamming_bound = Number[])

    for i in A
        push!(D, i)
    end

    CSV.write(joinpath(homedir(), "Desktop", "hamming_bound_integers.csv"), D)
end

# make_csv()

function bound_comparison()::Array{Array{Number, 1}}
    A = Array{AbstractFloat}[]
    upper_bound = 10
    increment_bound = 10
    stop_at = 10_000
    i = 0
    
    while true
        for q in 1:upper_bound, n in 1:upper_bound, d in 1:upper_bound
            if isprimepower(big(q))
                __arelessthan(upper_bound - increment_bound, q, n, d) && continue
                hb = hamming_bound(q, n, d)
                sb = singleton_bound(q, n, d)
                if ! isone(hb)
                    comparison = hb > sb ? 1 : (sb > hb ? -1 : 0)
                    push!(A, [q, n, d, hb, sb, comparison])
                    
                    i+=1
                    isequal(i, stop_at) && return A
                end
            end
        end
        upper_bound += increment_bound
    end
end

# A = bound_comparison()
# D = DataFrame(q = Number[], n = Number[], d = Number[], hamming_bound = Number[], singleton_bound = Number[], comparison = Number[])
#
# for i in A
#     push!(D, i)
# end
#
# CSV.write(joinpath(homedir(), "Desktop", "bound_comparison.csv"), D)
#
# theme(:solarized)
# variable_of_interest = :n
# plot = @df D scatter(
# 		variable_of_interest,
# 		:comparison,
#         # seriestype=:scatter,
# 		title = "",
# 		label = false,
# 		xaxis = "$(string(variable_of_interest))",
# 		yaxis = "Accuracy",
# 		fontfamily = font("Times"),
# 	)
#
# savefig(plot, joinpath(homedir(), "Desktop", "bound_comparison_by_$(string(variable_of_interest))_10000.pdf"))

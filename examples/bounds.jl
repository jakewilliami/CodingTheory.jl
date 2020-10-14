#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

include(joinpath(dirname(dirname(@__FILE__)), "src", "distance.jl"))
include(joinpath(dirname(dirname(@__FILE__)), "src", "utils.jl"))
include(joinpath(dirname(dirname(@__FILE__)), "src", "primes.jl"))

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

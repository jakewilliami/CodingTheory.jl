#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

using CSV, DataFrames, StatsPlots

include(joinpath(dirname(dirname(@__FILE__)), "src", "messages.jl"))
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

# You might want to use your code to search for parameters where q^n divided by sphere_packing_bound(q,n,d) is an integer. I think this is essentially how Golay started his search - it would be interesting to replicate some of that search.
function integer_search()::Array{Array{Number, 1}}
    A = []
    upper_bound = 10
    increment_bound = 10
    stop_at = 10_000
    
    while true
        for q in 1:upper_bound, n in 1:upper_bound, d in 1:upper_bound
			# skip configurations that have already been processed
			__arelessthan(upper_bound - increment_bound, q, n, d) && continue
            # note that we actually don't want to filter out non-linear codes
			# distance shouldn't be larger than the block length; filter trivial distances of one; we filter out combinations when all are equal; filter out codes that are generalised hamming codes; filter out even distances
			if d < n && ! isone(d) && ! __allequal(q, n, d) && ! ishammingbound(q, n, d) && ! iszero(mod(d, 2))
                hb = hamming_bound(q, n, d, no_round)
                if isinteger(hb) && ! isone(hb)
                    push!(A, [q, n, d, hb])
                    isequal(length(A), stop_at) && return A
                end
            end
        end
        upper_bound += increment_bound
    end
end

function make_integer_csv()
    A = integer_search()
    D = DataFrame(q = Number[], n = Number[], d = Number[], hamming_bound = Number[])

    for i in A
        push!(D, i)
    end

    CSV.write(joinpath(homedir(), "Desktop", "hamming_bound_integers.csv"), D)
    CSV.write(joinpath(dirname(dirname(@__FILE__)), "other", "hamming_bound_integers.csv"), D)
	println("Wrote data to $(joinpath(homedir(), "Desktop", "hamming_bound_integers.csv")).")
end

function bound_comparison()::Tuple{Int, Array{Array{Number, 1}}}
    A = Array{AbstractFloat}[]
    upper_bound = 10
    increment_bound = 10
    stop_at = 100_000
    
    while true
        for q in 1:upper_bound, n in 1:upper_bound, d in 1:upper_bound
			# skip configurations that have already been processed
			__arelessthan(upper_bound - increment_bound, q, n, d) && continue
			# remove q non prime powers
            if isprimepower(big(q))
                hb = hamming_bound(q, n, d)
                sb = singleton_bound(q, n, d)
				# filter out the more trivial bounds that are one (or if distance is one); filter out distances that are more than or equal to the block length, as they don't make sense in our context; filter out even distances
				# if d < n || 1 ∉ (hb, sb, d)
				if d < n && ! isone(d) && ! __allequal(q, n, d) && ! ishammingbound(q, n, d) && ! iszero(mod(d, 2))
                    comparison = hb > sb ? 1 : (sb > hb ? -1 : 0)
                    push!(A, [q, n, d, hb, sb, comparison])
                    isequal(length(A), stop_at) && return stop_at, A
                end
            end
        end
        upper_bound += increment_bound
    end
end

function plot_bound_comparison()
	number_of_bounds, A = bound_comparison()
	D = DataFrame(q = Number[], n = Number[], d = Number[], hamming_bound = Number[], singleton_bound = Number[], comparison = Number[])

	for i in A
	    push!(D, i)
	end

	CSV.write(joinpath(homedir(), "Desktop", "bound_comparison.csv"), D)
	println("Data written to $(joinpath(homedir(), "Desktop", "bound_comparison.csv"))")
	CSV.write(joinpath(dirname(dirname(@__FILE__)), "other", "bound_comparison.csv"), D)
	
	`Rscript "~"/projects/CodingTheory.jl/other/regression-tree.R`
	println("Regression tree saved at $(joinpath(dirname(dirname(@__FILE__)), "other", "Rplots.pdf"))")
end

# make_integer_csv constructs a csv file which looks for hamming bounds (given certain conditions) that are integers before rounding
# make_integer_csv()

# plot_bound_comparison constructs a dataframe which compares the hamming bound to the singleton bound and plots a regression tree using r
# plot_bound_comparison()


# If q is a prime power, and d = 3, and n is equal to q^((q^r-1)/(q-1)) for some positive integer r, then we know that the Hamming bound is integral, and we know that there is a perfect code with those parameters. So this choice of parameters is not terribly interesting.
#
# So I would filter as follows:

# This function tells us whether the Hamming bound we get for (n,d,q) is exactly the bound of a Hamming code. So if it returns true, then I would omit that triple (n,d,q) from the output.

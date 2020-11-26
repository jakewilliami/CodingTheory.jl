#!/usr/bin/env bash
    #=
    exec julia -tauto --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
	
include("plot_random_core.jl")

if isempty(ARGS) || length(ARGS) < 4
	throw(error("""
	Please specify the number of datapoints you want to calculate from the command line, along with q, n, and d.  For example
		./examples/plot_random.jl 6 7 3 1000 10000 # q, n, d, num_data_points, selection_set_size
	That is respectively: the size of the alphabet; the block length of the words; the distance of the code; the number of datapoints in the code; and the selection set size of each "random" sample.
	"""))
end

global q = parse(Int, ARGS[1])
global n = parse(Int, ARGS[2])
global d = parse(Int, ARGS[3])
global stop_at = parse(Int, ARGS[4])
global m = ifelse(length(ARGS) > 4, parse(Int, ARGS[5]), 10_000)

global upper_bound_adjustment = 10
global 🐖 = length(get_codewords_greedy(q, n, d))
global 🍖 = hamming_bound(q, n, d)
global 🐺 = singleton_bound(q, n, d)
graphing(q, n, d, stop_at, m = m)

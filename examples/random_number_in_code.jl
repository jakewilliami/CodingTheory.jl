#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

ENV["GKS_ENCODING"] = "utf-8"

using Plots, ProgressMeter, Formatting#, CSV, DataFrames

include(joinpath(dirname(dirname(@__FILE__)), "src", "CodingTheory.jl"))
using .CodingTheory

stop_at = parse(BigInt, ARGS[1])

q, n, d = 2, 7, 3

function graphing_points(stop_at::Integer)
    num_of_datapoints = stop_at
    data = []
	
	🐖 = length(get_codewords_greedy(q, n, d))
    
    @showprogress for _ in 1:num_of_datapoints
		push!(data, length(get_codewords_random(q, n, d)))
    end

    save_path = joinpath(dirname(@__DIR__), "other", "random_number_in_code,q=$(q),n=$n,d=$d,i=$num_of_datapoints,max.pdf")
    theme(:solarized)
	
	bins = minimum(data) + 0.5 : max(maximum(data), 🐖) + 0.5
        
    plot_points = histogram(
						data,
						fontfamily = font("Times"),
						xlabel = "|C|",
						ylabel = "Frequency of Times Found",
						label = "",
						title = "(q, n, d) = ($q, $n, $d): $(format(num_of_datapoints, commas = true)) Random Samples",
						xticks = minimum(data) : max(maximum(data), 🐖),
						bins = bins
						)
	
	🐖_idx = searchsortedlast(bins, 🐖) # <- the index of that value
	y_top = plot_points.series_list[1].plotattributes[:y][1:6:end][🐖_idx] # height of bar of choice
	
	annotate!(plot_points, 🐖, y_top + 1, "↓")

    savefig(plot_points, save_path)
    println("Plot saved at $save_path.")
end

graphing_points(stop_at)

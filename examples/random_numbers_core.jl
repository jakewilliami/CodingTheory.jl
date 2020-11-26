#!/usr/bin/env bash
    #=
    exec julia -tauto --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

println("Loading required libraries...\n")

using Base.Threads: @threads
using ProgressMeter: @showprogress, Progress, next!
using FLoops: @floop, ThreadedEx
using Plots, Formatting, CSV, DataFrames

plotly()

include(joinpath(dirname(dirname(@__FILE__)), "src", "CodingTheory.jl"))
using .CodingTheory

function obtain_data(q::Integer, n::Integer, d::Integer, stop_at::Integer; m::Integer=10_000)
	n < d && throw(error("You cannot have a distance greater than the block length."))
    data = []
	data_path = joinpath(dirname(@__DIR__), "other", "random_number_analysis", "random_number_in_code,q=$(q),n=$(n),d=$(d),i=$(stop_at),m=$(m).csv")
	
	## Obtain data
	p = Progress(stop_at, 1) # minimum update interval: 1 second
    @floop ThreadedEx() for _ in 1:stop_at
		push!(data, length(get_codewords_random(q, n, d, m = m)))
		next!(p) # increment progress bar
    end
	
	df = DataFrame(size_of_code = data)
	CSV.write(data_path, df)
	
	println("Wrote data to $(data_path)")
	
	return data
end

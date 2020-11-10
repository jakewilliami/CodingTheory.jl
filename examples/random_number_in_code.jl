#!/usr/bin/env bash
    #=
    exec julia -tauto --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

ENV["GKS_ENCODING"] = "utf-8"

using Base.Threads: @threads
using ProgressMeter: @showprogress, Progress, next!
using FLoops: @floop, ThreadedEx
using Plots, Formatting, CSV, DataFrames

plotly()

include(joinpath(dirname(dirname(@__FILE__)), "src", "CodingTheory.jl"))
using .CodingTheory

if isempty(ARGS)
	throw(error("Please specify the number of datapoints you want to calculate from the command line."))
end

stop_at = parse(BigInt, ARGS[1])

function add_plot_info(plt::AbstractPlot; extra::AbstractString=string(), counter::Integer=0)
    x = xlims(Plots.current())[1]
    y = ylims(Plots.current())[2]
    factor1 = y > 1 ? 2 : 1.6
    factor2 = y > 1 ? 4 : 2.5
    # plot!(annotations = (max(x * 2, x+100), (y - 10 * counter) / factor1,
    #        text(extra, :left, 13)))
	annotate!(plt, max(x * 2, x + 100), ((y - 10 * counter) / factor1), extra)
end

function obtain_data(q::Integer, n::Integer, d::Integer, stop_at::Integer; m::Integer=10_000)
	n < d && throw(error("You cannot have a distance greater than the block length."))
	upper_bound_adjustment = 10
    num_of_datapoints = stop_at
    data = []
	
	if @isdefined m
		batch_size = m
	else
		batch_size = 10_000 # the default value for m; a kwarg in `get_codewords_random`.
	end
	
	## Obtain data
	🐖 = length(get_codewords_greedy(q, n, d))
	🍖 = hamming_bound(q, n, d)
	🐺 = singleton_bound(q, n, d)
    
	p = Progress(num_of_datapoints, 1) # minimum update interval: 1 second
    @floop ThreadedEx() for _ in 1:num_of_datapoints
		push!(data, length(get_codewords_random(q, n, d, m = m)))
		next!(p) # increment progress bar
    end

    save_path = joinpath(dirname(@__DIR__), "other", "random_number_analysis", "random_number_in_code,q=$(q),n=$(n),d=$(d),i=$(num_of_datapoints),m=$(batch_size).pdf")
	data_path = joinpath(dirname(@__DIR__), "other", "random_number_analysis", "random_number_in_code,q=$(q),n=$(n),d=$(d),i=$(num_of_datapoints),m=$(batch_size).csv")
	
	D = DataFrame(size_of_code = Number[])

    for i in A
        push!(D, i)
    end
	
	CSV.write(data_path, D)
	
	println("Wrote data to $(data_path)")
	
	return data
end

function graphing(q::Integer, n::Integer, d::Integer, stop_at::Integer; m::Integer=10_000)
	data = obtain_data(q, n, d, stop_at, m = m)
	## Plot points
	theme(:solarized)
	
	x_max = ifelse(abs(🐖 - 🍖) ≤ upper_bound_adjustment && abs(🐖 - 🐺) ≤ upper_bound_adjustment, maximum([🐖, 🍖, 🐺]), 🐖)
	# x_max = ifelse(abs(🐖 - 🍖) ≤ upper_bound_adjustment && abs(🐖 - 🐺) ≤ upper_bound_adjustment, maximum([🐖, 🍖, 🐺]), ifelse(abs(🐖 - 🍖) ≤ upper_bound_adjustment, max(🐖, 🍖), ifelse(abs(🐖 - 🐺) ≤ upper_bound_adjustment, max(🐖, 🍖), 🐖)))
	bins = minimum(data) - 0.5 : max(maximum(data), x_max) + 0.5
    
    plt = histogram(
					data,
					fontfamily = font("Times"),
					xlabel = "|C|",
					ylabel = "Frequency of Times Found",
					label = "",
					title = "(q, n, d) = ($q, $n, $d): $(format(num_of_datapoints, commas = true)) Random Samples, with Selection Size $(format(batch_size, commas = true))",
					xticks = minimum(data) : max(maximum(data), x_max),
					bins = bins
					)
	
	## Annotate plot
	annotation_counter = 0
	if isempty(data)
		add_plot_info(plt, extra="🐖 = $(🐖)", counter=annotation_counter); annotation_counter += 1
		add_plot_info(plt, extra="🍖 = $(🍖)", counter=annotation_counter); annotation_counter += 1
		add_plot_info(plt, extra="🐺 = $(🐺)", counter=annotation_counter); annotation_counter += 1
		savefig(plt, save_path)
	    println("Plot saved at $save_path.")
		return nothing
	end

	🐖_idx = searchsortedlast(bins, 🐖) # <- the index of that value
	🐖_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🐖_idx] # height of bar of choice
	annotate!(plt, 🐖, 🐖_y_top + 1, "🐖")
	
	🍖_y, 🐺_y = 0, 0
	🍖_annotated, 🐺_annotated = false, false
	
	if abs(🐖 - 🍖) ≤ upper_bound_adjustment && abs(🐖 - 🐺) ≤ upper_bound_adjustment # the singleton AND hamming bounds can fit on the plot
		🍖_idx = searchsortedlast(bins, 🍖) # <- the index of that value
		🍖_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🍖_idx] # height of bar of choice
		🐺_idx = searchsortedlast(bins, 🐺) # <- the index of that value
		🐺_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🐺_idx] # height of bar of choice
	
		# adjust emoticon positioning if overlapping values.
		if 🐖 == 🍖
			🍖_y = 🍖_y_top + 10
		end
	
		if 🐖 == 🐺
			if 🐖 == 🍖 # then all of 🍖, 🐖, and 🐺 are the same
				🍖_y = 🍖_y_top + 10
				🐺_y = 🐺_y_top + 20
			else
				🐺_y = 🐺_y_top + 10
			end
		end
	
		annotate!(plt, 🍖, 🍖_y, "🍖")
		annotate!(plt, 🐺, 🐺_y, "🐺")
			
		🍖_annotated = true
		🐺_annotated = true
	end
	
	if abs(🐖 - 🍖) ≤ upper_bound_adjustment && ! 🍖_annotated# only the hamming bound fits on the plot
		🍖_idx = searchsortedlast(bins, 🍖) # <- the index of that value
		🍖_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🍖_idx] # height of bar of choice
	
		# adjust emoticon positioning if overlapping values.
		if 🐖 == 🍖
			🍖_y = 🍖_y_top + 10
		end
	
		annotate!(plt, 🍖, 🍖_y, "🍖")
	
		🍖_annotated = true
	end
	
	if abs(🐖 - 🐺) ≤ upper_bound_adjustment && ! 🐺_annotated # only the singleton bound fits on the plot
		🐺_idx = searchsortedlast(bins, 🐺) # <- the index of that value
		🐺_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🐺_idx] # height of bar of choice
	
		# adjust emoticon positioning if overlapping values.
		if 🐖 == 🐺
			🐺_y = 🐺_y_top + 10
		end
	
		annotate!(plt, 🐺, 🐺_y, "🐺")
	
		🐺_annotated = true
	end
	
	if ! 🍖_annotated
		annotate!(plt, maximum(filter(!isnan, plt.series_list[1].plotattributes[:x])), maximum(filter(!isnan, plt.series_list[1].plotattributes[:y])) + 0, "🍖 = $(🍖)")
		# add_plot_info(plt, extra="🍖 = $(🍖)", counter=annotation_counter); annotation_counter += 1
	end
	
	if ! 🐺_annotated
		annotate!(plt, maximum(filter(!isnan, plt.series_list[1].plotattributes[:x])), maximum(filter(!isnan, plt.series_list[1].plotattributes[:y])) - 10, "🐺 = $(🐺)")
		# add_plot_info(plt, extra="🐺 = $(🐺)", counter=annotation_counter); annotation_counter += 1
	end

    savefig(plt, save_path)
    println("Plot saved at $save_path.")
	return nothing
end

# for i in 3:7
global q, n, d = 2, 2, 2
graphing(q, n, d, stop_at, m = 100)
# end

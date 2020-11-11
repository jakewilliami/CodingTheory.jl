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

function graphing(q::Integer, n::Integer, d::Integer, data::Union{AbstractArray, DataFrame}, stop_at::Integer; m::Integer=10_000)
	save_path = joinpath(dirname(@__DIR__), "other", "random_number_analysis", "random_number_in_code,q=$(q),n=$(n),d=$(d),i=$(stop_at),m=$(m).pdf")
	
	if data isa DataFrame
		data = Matrix(data)
	end
	
	🐖_obtained = ifelse(🐖 ∈ data, true, false)
	
	## Plot points
	theme(:solarized)
	
	x_max = 0
	if 🍖 - maximum(data) ∈ [1:upper_bound_adjustment...]
		x_max = 🍖
	end
	if 🐺 - maximum(data) ∈ [1:upper_bound_adjustment...]
		x_max = 🐺
	end
	bin_adjustment = ifelse(🐖_obtained, 0.5, minimum(abs.(extrema(vcat(data, x_max)) .- 🐖)) + 0.5) # ensure the bins encompass the greedy
	bins = minimum(data) - bin_adjustment : max(maximum(data), x_max) + bin_adjustment
    
    plt = histogram(
					data,
					fontfamily = font("Times"),
					xlabel = "|C|",
					ylabel = "Frequency of Times Found",
					label = "",
					title = "(q, n, d) = ($q, $n, $d): $(format(stop_at, commas = true)) Random Samples, with Selection Size $(format(m, commas = true))",
					# xticks = minimum(data) : max(maximum(data), x_max),
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
	annotate!(plt, 🐖, 🐖_y_top, "🐖")
	
	annotation_factor = 10
	x_annotation_factor = abs(bins[end] - bins[1]) * 0.05
	y_annotation_factor = abs(bins[end] - bins[1]) * 1.2
	🍖_y, 🐺_y = 0, 0
	🍖_annotated, 🐺_annotated = false, false
	
	if abs(x_max - 🍖) ≤ upper_bound_adjustment && abs(x_max - 🐺) ≤ upper_bound_adjustment # the singleton AND hamming bounds can fit on the plot
		🍖_idx = searchsortedlast(bins, 🍖) # <- the index of that value
		🍖_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🍖_idx] # height of bar of choice
		🐺_idx = searchsortedlast(bins, 🐺) # <- the index of that value
		🐺_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🐺_idx] # height of bar of choice
	
		# adjust emoticon positioning if overlapping values.
		if 🐖 == 🍖
			🍖_y = 🍖_y_top + y_annotation_factor
		end
	
		if 🐖 == 🐺
			if 🐖 == 🍖 # then all of 🍖, 🐖, and 🐺 are the same
				🍖_y = 🍖_y_top + y_annotation_factor
				🐺_y = 🐺_y_top + 2*y_annotation_factor
			else
				🐺_y = 🐺_y_top + y_annotation_factor
			end
		end
	
		annotate!(plt, 🍖, 🍖_y, "🍖")
		annotate!(plt, 🐺, 🐺_y, "🐺")
			
		🍖_annotated = true
		🐺_annotated = true
	end
	
	if abs(x_max - 🍖) ≤ upper_bound_adjustment && ! 🍖_annotated # only the hamming bound fits on the plot
		🍖_idx = searchsortedlast(bins, 🍖) # <- the index of that value
		🍖_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🍖_idx] # height of bar of choice
	
		# adjust emoticon positioning if overlapping values.
		if 🐖 == 🍖
			🍖_y = 🍖_y_top + y_annotation_factor
		end
	
		annotate!(plt, 🍖, 🍖_y, "🍖")
	
		🍖_annotated = true
	end
	println(🐺)
	println(abs(x_max - 🐺) ≤ upper_bound_adjustment)
	if abs(x_max - 🐺) ≤ upper_bound_adjustment && ! 🐺_annotated # only the singleton bound fits on the plot
		🐺_idx = searchsortedlast(bins, 🐺) # <- the index of that value
		🐺_y_top = plt.series_list[1].plotattributes[:y][1:6:end][🐺_idx] # height of bar of choice
	
		# adjust emoticon positioning if overlapping values.
		if 🐖 == 🐺
			🐺_y = 🐺_y_top + y_annotation_factor
		end
	
		annotate!(plt, 🐺, 🐺_y, "🐺")
	
		🐺_annotated = true
	end
	
	if ! 🍖_annotated
		annotate!(plt, maximum(filter(!isnan, plt.series_list[1].plotattributes[:x])) - x_annotation_factor, maximum(filter(!isnan, plt.series_list[1].plotattributes[:y])) + 0, "🍖 = $(🍖)")
		# add_plot_info(plt, extra="🍖 = $(🍖)", counter=annotation_counter); annotation_counter += 1
	end
	
	if ! 🐺_annotated
		annotate!(plt, maximum(filter(!isnan, plt.series_list[1].plotattributes[:x])) - x_annotation_factor, maximum(filter(!isnan, plt.series_list[1].plotattributes[:y])) - y_annotation_factor, "🐺 = $(🐺)")
		# add_plot_info(plt, extra="🐺 = $(🐺)", counter=annotation_counter); annotation_counter += 1
	end

    savefig(plt, save_path)
    println("Plot saved at $save_path.")
	return nothing
end

graphing(q::Integer, n::Integer, d::Integer, stop_at::Integer; m::Integer=10_000) =
	graphing(q, n, d, obtain_data(q, n, d, stop_at, m = m), stop_at, m = m)

# if isempty(ARGS) || length(ARGS) < 4
# 	throw(error("""
# 	Please specify the number of datapoints you want to calculate from the command line, along with q, n, and d.  For example
# 		./examples/random_number_in_code.jl 6 7 3 1000 10000
# 	That is respectively: the size of the alphabet; the block length of the words; the distance of the code; the number of datapoints in the code; and the selection set size of each "random" sample.
# 	"""))
# end

# global q = parse(Int, ARGS[1])
# global n = parse(Int, ARGS[2])
# global d = parse(Int, ARGS[3])
# global stop_at = parse(Int, ARGS[4])
# global m = ifelse(length(ARGS) > 4, parse(Int, ARGS[5]), 10_000)
#
# global upper_bound_adjustment = 10
# global 🐖 = length(get_codewords_greedy(q, n, d))
# global 🍖 = hamming_bound(q, n, d)
# global 🐺 = singleton_bound(q, n, d)




# (6, 3, 2),
# (10, 3, 2),
# for (q, n, d) in [(6, 4, 2), (12, 3, 2), (21, 4, 3), (14, 3, 2), (15, 3, 2), (18, 3, 2), (6, 5, 2), (6, 7, 3)]
# 	println("Processing (q, n, d) = ($(q), $(n), $(d))")
# 	global upper_bound_adjustment = 10
# 	global 🐖 = length(get_codewords_greedy(q, n, d))
# 	global 🍖 = hamming_bound(q, n, d)
# 	global 🐺 = singleton_bound(q, n, d)
#
# 	graphing(q, n, d, 1000, m = 10_000)
# end



# graphing(q, n, d, stop_at, m = m)

# reading from pre-obtained data
q, n, d = 6, 3, 2
m = 10_000
stop_at = 1000
global upper_bound_adjustment = 10
global 🐖 = length(get_codewords_greedy(q, n, d))
global 🍖 = hamming_bound(q, n, d)
global 🐺 = singleton_bound(q, n, d)
# data = CSV.read("/Users/jakeireland/projects/CodingTheory.jl/other/random_number_analysis/random_number_in_code,q=$(q),n=$(n),d=$(d),i=$(stop_at),m=$(m).csv")
# data = CSV.read("/Users/jakeireland/projects/CodingTheory.jl/other/random_number_analysis/mildly_interesting_ones/[2]-10-3-2/random_number_in_code,q=10,n=3,d=2,i=1000,m=10000.csv")
data = CSV.read("/Users/jakeireland/projects/CodingTheory.jl/other/random_number_analysis/mildly_interesting_ones/[1]-6-2-3/random_number_in_code,q=6,n=3,d=2,i=1000,m=10000.csv")
graphing(q, n, d, data, stop_at, m = m)

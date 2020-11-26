#!/usr/bin/env bash
    #=
    exec julia -tauto --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
ENV["GKS_ENCODING"] = "utf-8"

include("random_numbers_core.jl")

function add_plot_info(plt::AbstractPlot; extra::AbstractString=string(), counter::Integer=0)
    x = xlims(Plots.current())[1]
    y = ylims(Plots.current())[2]
    factor1 = y > 1 ? 2 : 1.6
    factor2 = y > 1 ? 4 : 2.5
    # plot!(annotations = (max(x * 2, x+100), (y - 10 * counter) / factor1,
    #        text(extra, :left, 13)))
	annotate!(plt, max(x * 2, x + 100), ((y - 10 * counter) / factor1), extra)
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
	if 🍖 - maximum(data) ∈ [0:upper_bound_adjustment...]
		x_max = 🍖
	end
	if 🐺 - maximum(data) ∈ [0:upper_bound_adjustment...]
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

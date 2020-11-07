#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

using StatsPlots, ProgressMeter, CSV, DataFrames

include(joinpath(dirname(dirname(@__FILE__)), "src", "CodingTheory.jl"))
using .CodingTheory

stop_at = parse(BigInt, ARGS[1])

function graphing_time(stop_at::Integer)
    stop_q_at = 8
    stop_n_at = 8
    num_of_datapoints = stop_at
    data = Matrix{Float64}(undef, num_of_datapoints, 3) # needs 3 columns: q, n, and time
    
    @showprogress for i in 1:num_of_datapoints
        q = rand(1:stop_q_at)
        n = rand(1:stop_n_at)
        time_took = @elapsed get_all_words(q, n)
        data[i, :] .= [q, n, time_took]
    end

    save_path = joinpath(dirname(@__DIR__), "other", "get_codewords_time_analysis,n=$stop_n_at,m=$stop_q_at,i=$num_of_datapoints.pdf")
    theme(:solarized)
        
    plot_points = scatter(data[:,1:2], data[:,3], fontfamily = font("Times"), ylabel = "Time elapsed during calculation [seconds]", clabel = "Value of q or n", label = ["q" "n"], smooth = true)#, xlims = (0, stop_n_at)) # smooth = true
    
    savefig(plot_points, save_path)
    println("Plot saved at $save_path.")
    
    D = DataFrame(q = Number[], n = Number[], time = Number[])

    for i in eachrow(data)
        push!(D, i)
    end
	
	data_file_desktop = joinpath(homedir(), "Desktop", "get_codewords_time_analysis,n=$stop_n_at,m=$stop_q_at,i=$num_of_datapoints.csv")
	data_file_other = joinpath(dirname(dirname(@__FILE__)), "other", "get_codewords_time_analysis,n=$stop_n_at,m=$stop_q_at,i=$num_of_datapoints.csv")

    CSV.write(data_file_desktop, D)
    CSV.write(data_file_other, D)
	println("Wrote data to $(data_file_other).")
end

graphing_time(stop_at)

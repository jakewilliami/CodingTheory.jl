#=
Example usage:
```
	./examples/plot_random_from_file.jl /Users/jakeireland/Desktop/random_number_in_code,q=21,n=4,d=3,i=1000,m=10000.csv 6 7 3 # dataset, q, n, d
```
=#

include("plot_random_core.jl")

# reading from pre-obtained data
if isempty(ARGS) || length(ARGS) < 5
    throw(error("""
    Please specify the dataset from the command line, along with q, n, and d.  For example
    ```
    	./examples/plot_random_from_file.jl /Users/jakeireland/Desktop/random_number_in_code,q=21,n=4,d=3,i=1000,m=10000.csv 6 7 3 10000# dataset, q, n, d, selection_set_size
    ```
    That is respectively: the size of the alphabet; the block length of the words; the distance of the code; the number of datapoints in the code; and the selection set size of each "random" sample.
    """,),)
end

global data = ARGS[1] # e.g. "/Users/jakeireland/Desktop/random_number_in_code,q=21,n=4,d=3,i=1000,m=10000.csv"
global q = parse(Int, ARGS[2])
global n = parse(Int, ARGS[3])
global d = parse(Int, ARGS[4])
global m = parse(Int, ARGS[5])

global upper_bound_adjustment = 10
global 🐖 = length(get_codewords_greedy(q, n, d))
global 🍖 = hamming_bound(q, n, d)
global 🐺 = singleton_bound(q, n, d)
# data = CSV.read("/Users/jakeireland/projects/CodingTheory.jl/other/random_number_analysis/random_number_in_code,q=$(q),n=$(n),d=$(d),i=$(stop_at),m=$(m).csv")
# data = CSV.read("/Users/jakeireland/projects/CodingTheory.jl/other/random_number_analysis/mildly_interesting_ones/[2]-10-3-2/random_number_in_code,q=10,n=3,d=2,i=1000,m=10000.csv")
data = CSV.read(data)
graphing(q, n, d, data, nrow(data); m = m)

using DataFrames, CSV

# use bounds.jl to get the hamming bounds
datafile = "/Users/jakeireland/Desktop/hamming_bound_integers_10000.csv"

df = sort(DataFrame(CSV.read(datafile)), :hamming_bound)
df[!,:] = round.(BigInt, df[!,:])

println(df[1:20, :])


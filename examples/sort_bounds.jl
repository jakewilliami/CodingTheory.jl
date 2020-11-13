#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

include(joinpath(dirname(@__DIR__), "src", "CodingTheory.jl"))
using .CodingTheory

using DataFrames, CSV
using DataFramesMeta: @where
# using Lazy: @>

const ruled_out = Tuple[
        # (3, 11, 5), # n words ≈ 186 ≠ 729
        # (21, 4, 3), # n words ≈ 359 ≠ 2401
        # (6, 7, 3), # n words ≈ 2995 ≠ 7776
        # (2, 23, 7), # this is one of the Golay codes, silly me
    ]

function load_df(pathname::AbstractString)
    # use bounds.jl to get the hamming bounds
    df = DataFrame(CSV.read(pathname))
    insert!(df, ncol(df) + 1, CodingTheory.sizeof_perfect_code.(df.q, df.n, df.d), :size_of_perfect_code)
    df = round.(BigInt, df)

    return df
end

function filter_df(
    df::DataFrame;
    d_not_one::Bool=true,
    not_all_equal::Bool=true,
    not_hamming_perfect::Bool=true,
    d_not_even::Bool=true,
    q_prime_power::Bool=false,
    q_not_prime_power::Bool=true,
    hamming_bound_not_one::Bool=true,
    hamming_bound_is_lower_bound::Bool=true
    )

    if d_not_one
        df = @where(df, :d .≠ 1)
    end
    
    if not_all_equal
        df = @where(df, .! CodingTheory.allequal.(:q, :n, :d))
    end
    
    if not_hamming_perfect
        df = @where(df, .! ishammingperfect.(:q, :n, :d))
    end
    
    if d_not_even
        df = @where(df, .! iseven.(:d))
    end
    
    if q_prime_power
        df = @where(df, CodingTheory.isprimepower.(:q))
    end
    
    if q_not_prime_power
        df = @where(df, .! CodingTheory.isprimepower.(:q))
    end
    
    if hamming_bound_not_one
        df = @where(df, :hamming_bound .≠ 1)
    end

    if hamming_bound_is_lower_bound
        df = @where(df, isequal.(:hamming_bound, :smallest_bound))
    end
    
    # Filter out those (q, n, d) configurations already ruled out.
    # df = @where(df, tuple.(:q, :n, :d) .∉ Ref(ruled_out))
    
    return df
end

function view_df(
    df::DataFrame;
    nrows::Integer=10,
    d_not_one::Bool=true,
    not_all_equal::Bool=true,
    not_hamming_perfect::Bool=true,
    d_not_even::Bool=true,
    q_prime_power::Bool=false,
    q_not_prime_power::Bool=true,
    hamming_bound_not_one::Bool=true,
    hamming_bound_is_lower_bound::Bool=true
    )
    df = filter_df(
        df;
        d_not_one=d_not_one,
        not_all_equal=not_all_equal,
        not_hamming_perfect=not_hamming_perfect,
        d_not_even=d_not_even,
        q_prime_power=q_prime_power,
        q_not_prime_power=q_not_prime_power,
        hamming_bound_not_one=hamming_bound_not_one,
        hamming_bound_is_lower_bound=hamming_bound_is_lower_bound
        )
    
    # sort last!
    df = sort(df, :hamming_bound)
    df = sort(df, :size_of_perfect_code)
    println(size(df))
    return select(df, Not([:smallest_bound]))[1:nrows, :]
end





##############################################################################################################

# datafile = "/Users/jakeireland/Desktop/bounds/hamming_bound_integers_10000_incl_sb.csv" # obtain using bounds.jl
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_excl_sb.csv" # obtain using bounds.jl
# datafile = "/Users/jakeireland/Desktop/bounds/hamming_bound_integers_2000000_og.csv" # obtain using bounds.jl
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_slightly_looser.csv"
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_very_loose.csv"
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_very_loose_q_neq_1.csv"
datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_100000_very_loose_q_neq_1.csv" # <=====

df = view_df(
    load_df(datafile);
    nrows=10,
    d_not_one=true, # distance equal to one is trivially the list of all words with q and n
    not_all_equal=false,
    not_hamming_perfect=true,
    d_not_even=true, # cannot get a perfect code if d is even, because of how hamming bound is calculated
    q_prime_power=true, # inverse of the below
    q_not_prime_power=false, # we are interested in alphabet size not prime power
    hamming_bound_not_one=false,
    hamming_bound_is_lower_bound=true # cannnot obtain perfect code if singleton bound is less than hamming bound
    )

println(df)

df2 = filter_df(
    load_df(datafile);
    d_not_one=true,
    not_all_equal=false,
    not_hamming_perfect=true,
    d_not_even=true,
    q_prime_power=true,
    q_not_prime_power=false,
    hamming_bound_not_one=false,
    hamming_bound_is_lower_bound=true
    )
    
df2 = select(df2, Not([:smallest_bound]))

D = DataFrame(q = Number[], n = Number[], d = Number[], hamming_bound = Number[], size_of_perfect_code = Number[])

for row in eachrow(df2)
    if isinteger(log(row.q, hamming_bound(row.q, row.n, row.d)))
        push!(D, row)
    end
end

println(size(D))
println(sort(D, :size_of_perfect_code))

## original code:
# df = round.(BigInt, DataFrame(CSV.read(datafile)))
# df = @where(df, .! isprimepower.(:q))
# df = sort(df, :hamming_bound)[1:10, :]
# println(df)

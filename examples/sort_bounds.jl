using CodingTheory
using DataFrames, CSV
using DataFramesMeta: @where
using Formatting
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
    println("The complete list has $(format(nrow(df), commas = true)) observations.")
    insert!(
        df,
        ncol(df) + 1,
        CodingTheory.sizeof_perfect_code.(df.q, df.n, df.d),
        :size_of_perfect_code,
    )
    df = round.(BigInt, df)

    return df
end

function filter_df(
    df::DataFrame;
    q_not_one::Bool = true,
    d_not_one::Bool = true,
    d_not_two::Bool = true,
    d_not_n::Bool = true,
    d_leq_n::Bool = true,
    not_hamming_perfect::Bool = true,
    d_not_even::Bool = true,
    q_prime_power::Bool = false,
    q_not_prime_power::Bool = true,
    hamming_bound_isinteger::Bool = true,
    hamming_bound_is_lower_bound::Bool = true,
)
    if q_not_one
        df = @where(df, :q .≠ 1)
    end

    if d_not_one
        df = @where(df, :d .≠ 1)
    end

    if d_not_two
        df = @where(df, :d .≠ 2)
    end

    if d_not_n
        df = @where(df, :d .≠ :n)
    end

    if d_leq_n
        df = @where(df, :d .≤ :n)
    end

    if not_hamming_perfect
        df = @where(df, .!ishammingperfect.(:q, :n, :d))
    end

    if d_not_even
        df = @where(df, .!iseven.(:d))
    end

    if q_prime_power
        df = @where(df, CodingTheory.isprimepower.(:q))
    end

    if q_not_prime_power
        df = @where(df, .!CodingTheory.isprimepower.(:q))
    end

    if hamming_bound_isinteger
        df = @where(df, isinteger.(:hamming_bound))
    end

    if hamming_bound_is_lower_bound
        df = @where(df, isequal.(:hamming_bound, :smallest_bound))
        # df = @where(df, isequal.(:hamming_bound, min.(:hamming_bound, :singleton_bound)))
    end

    # Filter out those (q, n, d) configurations already ruled out.
    # df = @where(df, tuple.(:q, :n, :d) .∉ Ref(ruled_out))

    return df
end

function view_df(
    df::DataFrame;
    nrows::Integer = 10,
    q_not_one::Bool = true,
    d_not_one::Bool = true,
    d_not_two::Bool = true,
    d_not_n::Bool = true,
    d_leq_n::Bool = true,
    not_hamming_perfect::Bool = true,
    d_not_even::Bool = true,
    q_prime_power::Bool = false,
    q_not_prime_power::Bool = true,
    hamming_bound_isinteger::Bool = true,
    hamming_bound_is_lower_bound::Bool = true,
)
    df = filter_df(
        df;
        q_not_one = q_not_one,
        d_not_one = d_not_one,
        d_not_two = d_not_two,
        d_not_n = d_not_n,
        d_leq_n = d_leq_n,
        not_hamming_perfect = not_hamming_perfect,
        d_not_even = d_not_even,
        q_prime_power = q_prime_power,
        q_not_prime_power = q_not_prime_power,
        hamming_bound_isinteger = hamming_bound_isinteger,
        hamming_bound_is_lower_bound = hamming_bound_is_lower_bound,
    )

    # sort last!
    df = sort(df, :hamming_bound)
    # df = sort(df, :size_of_perfect_code)

    println("The filtered list has $(format(nrow(df), commas = true)) observations matching the filtered criteria.  Now showing the 10 smallest values.",)
    # return select(df, Not([:smallest_bound]))[1:nrows, :]
    return df[1:nrows, :]
    # return select(df, Not([:singleton_bound]))[1:nrows, :]
end

##############################################################################################################

datafile = "/Users/jakeireland/projects/CodingTheory.jl/other/research/bounds/bound_integers_100000_very_loose_q_neq_1.csv" # <=====
# datafile = "/Users/jakeireland/projects/CodingTheory.jl/other/bound_integers_1000000.csv"

# df = load_df("/Users/jakeireland/projects/CodingTheory.jl/other/research/bounds/bound_integers_100000_very_loose_q_neq_1.csv")
# println(df[end, :])

# df = load_df("/Users/jakeireland/projects/CodingTheory.jl/other/bound_integers_1000000.csv")

# println(df[end, :])

df = view_df(
    load_df(datafile);
    nrows = 10,
    q_not_one = true,
    d_not_one = true, # distance equal to one is trivially the list of all words with q and n
    d_not_two = true, # distance equal to two is trivial as the words' coordinates sum to zero modulo q
    d_not_n = true, # distance equal to two is trivially the repetition code
    d_leq_n = true,
    not_hamming_perfect = true,
    d_not_even = false, # cannot get a perfect code if d is even, because of how hamming bound is calculated
    q_prime_power = false, # inverse of the below
    q_not_prime_power = true, # we are interested in alphabet size not prime power
    hamming_bound_isinteger = true,
    hamming_bound_is_lower_bound = true, # cannnot obtain perfect code if singleton bound is less than hamming bound
)

println(df)

# df2 = filter_df(
#     load_df(datafile);
#     nrows=10,
#     q_not_one=true,
#     d_not_one=true, # distance equal to one is trivially the list of all words with q and n
#     d_not_two=true, # distance equal to two is trivial as the words' coordinates sum to zero modulo q
#     d_not_n=true, # distance equal to two is trivially the repetition code
#     d_leq_n=true,
#     not_hamming_perfect=true,
#     d_not_even=false, # cannot get a perfect code if d is even, because of how hamming bound is calculated
#     q_prime_power=false, # inverse of the below
#     q_not_prime_power=false, # we are interested in alphabet size not prime power
#     hamming_bound_isinteger=true,
#     hamming_bound_is_lower_bound=false # cannnot obtain perfect code if singleton bound is less than hamming bound
#     )
#
# df2 = select(df2, Not([:smallest_bound]))
#
# D = DataFrame(q = Number[], n = Number[], d = Number[], hamming_bound = Number[], size_of_perfect_code = Number[])
#
# for row in eachrow(df2)
#     if isinteger(log(row.q, hamming_bound(row.q, row.n, row.d)))
#         push!(D, row)
#     end
# end
#
# println(size(D))
# println(sort(D, :size_of_perfect_code))

## original code:
# df = round.(BigInt, DataFrame(CSV.read(datafile)))
# df = @where(df, .! isprimepower.(:q))
# df = sort(df, :hamming_bound)[1:10, :]
# println(df)

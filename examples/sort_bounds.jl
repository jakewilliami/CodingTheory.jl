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

#=
\section{Theorem} \emph{Let $q$ be a prime power and let $\Sigma$ be an alphabet of size $q$.  Let $C\subseteq\Sigma^n$ be a code with block length $n$, distance $d$, and $M$ codewords.  If $C$ is a non-trivial perfect code (i.e., $M$ is equal to the Hamming bound), then either:}\begin{enumerate}\item\emph{$n=\frac{q^r-1}{q-1}$ for some positive integer $r$, $d=3$, and $M=q^{\sfrac{q^r-1}{q-1}-r}$}\item\emph{$q=2$, $n=23$, $d=7$, and $M=2^{12}$, or}\item\emph{$q=3$, $n=11$, $d=5$, and $M=2^6$.}\end{enumerate}  This theorem shows that the only possible parameters for non-trivial perfect codes are the same as the parameters for Hamming and Golay codes, \emph{as long as the size of the alphabet is a prime power}.  This obviously leads to the following problem, which is still unanswered.  \section{Conjecture}\emph{Let $q$ be a positive integer that is not a prime power.  There does not exist a non-trivial perfect code with an alphabet of size $q$.}
=#

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
        q_not_prime_power=q_not_prime_power,
        hamming_bound_not_one=hamming_bound_not_one,
        hamming_bound_is_lower_bound=hamming_bound_is_lower_bound
        )
    
    # sort last!
    df = sort(df, :hamming_bound)
    df = sort(df, :size_of_perfect_code)

    return select(df, Not([:smallest_bound]))[1:nrows, :]
end





##############################################################################################################

# datafile = "/Users/jakeireland/Desktop/bounds/hamming_bound_integers_10000_incl_sb.csv" # obtain using bounds.jl
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_excl_sb.csv" # obtain using bounds.jl
# datafile = "/Users/jakeireland/Desktop/bounds/hamming_bound_integers_2000000_og.csv" # obtain using bounds.jl
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_slightly_looser.csv"
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_very_loose.csv"
# datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_10000_very_loose_q_neq_1.csv"
datafile = "/Users/jakeireland/Desktop/bounds/bound_integers_100000_very_loose_q_neq_1.csv"

df = view_df(
    load_df(datafile);
    nrows=10,
    d_not_one=true,
    not_all_equal=false,
    not_hamming_perfect=false,
    d_not_even=false,
    q_not_prime_power=true,
    hamming_bound_not_one=false,
    hamming_bound_is_lower_bound=true
    )

println(df)


## original code:
# df = round.(BigInt, DataFrame(CSV.read(datafile)))
# df = @where(df, .! isprimepower.(:q))
# df = sort(df, :hamming_bound)[1:10, :]
# println(df)

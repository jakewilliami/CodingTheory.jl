using DataFrames, CSV
using DataFramesMeta: @where
# using Lazy: @>

#=
\section{Theorem} \emph{Let $q$ be a prime power and let $\Sigma$ be an alphabet of size $q$.  Let $C\subseteq\Sigma^n$ be a code with block length $n$, distance $d$, and $M$ codewords.  If $C$ is a non-trivial perfect code (i.e., $M$ is equal to the Hamming bound), then either:}\begin{enumerate}\item\emph{$n=\frac{q^r-1}{q-1}$ for some positive integer $r$, $d=3$, and $M=q^{\sfrac{q^r-1}{q-1}-r}$}\item\emph{$q=2$, $n=23$, $d=7$, and $M=2^{12}$, or}\item\emph{$q=3$, $n=11$, $d=5$, and $M=2^6$.}\end{enumerate}  This theorem shows that the only possible parameters for non-trivial perfect codes are the same as the parameters for Hamming and Golay codes, \emph{as long as the size of the alphabet is a prime power}.  This obviously leads to the following problem, which is still unanswered.  \section{Conjecture}\emph{Let $q$ be a positive integer that is not a prime power.  There does not exist a non-trivial perfect code with an alphabet of size $q$.}
=#

const ruled_out = Tuple[
        (3, 11, 5), # n words ≈ 186 ≠ 729
        (21, 4, 3), # n words ≈ 359 ≠ 2401
        (6, 7, 3), # n words ≈ 2995 ≠ 7776
        # (2, 23, 7), # this is one of the Golay codes, silly me
    ]



# use bounds.jl to get the hamming bounds
datafile = "/Users/jakeireland/Desktop/hamming_bound_integers_10000.csv" # obtain using bounds.jl

# Sort dataframe by hamming bound
df = round.(BigInt, sort(DataFrame(CSV.read(datafile)), :hamming_bound))

# println(df[1:11, :])

# Filter out those (q, n, d) configurations already ruled out.
df = @where(df, tuple.(:q, :n, :d) .∉ Ref(ruled_out))

println(df[1:11, :])

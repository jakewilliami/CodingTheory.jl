#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
module CodingTheory

ENV["HECKE_PRINT_BANNER"] = "false"
# ENV["NEMO_PRINT_BANNER"] = "false"

# using Nemo: isprime, factor, fmpz # for prime power function
using Hecke: ispower, isprime_power, isprime # for perfect power function
# using Primes: isprime, primes
using LinearAlgebra: I
using FLoops: @floop, ThreadedEx
using Polynomials, StaticArrays

# from utils.jl
export sizeof_perfect_code, sizeof_all_words, has_identity, has_identity_on_left

include("utils.jl")

# Abstract types
export FinitePolynomial, AbstractCode, NonStaticAbstractWord, AbstractWord, Word,
        Codewords, Alphabet, Messages, CodeUniverse, CodeUniverseIterator,
        UniverseParameters
export no_round, getindex, setindex!, firstindex, lastindex, size, length, rand,
        gensym, genalphabet, eltype, isword, isabstractword

# RREF
export rref, rref!

# Bounds, Messages, Distance, and Primes
export rate, sphere_covering_bound, sphere_packing_bound, hamming_bound,
        singleton_bound, gilbert_varshamov_bound, elias_bassalygo_bound,
        plotkin_bound, construct_ham_matrix, isperfect, ishammingperfect,
        isgolayperfect
export get_codewords_greedy, get_codewords_random, get_all_words, get_codewords
export hamming_distance, hamming_ball, code_distance, t_error_detecting,
        t_error_correcting, find_error_detection_max, find_error_correction_max
export isprimepower

# Algebra
export FinitePolynomial, list_polys, multiplication_table, list_span, islinear,
        isirreducible, normal_form!, normal_form, equivalent_code!, equivalent_code,
        generator!, generator, parity_check, syndrome, isincode

# Powers
export isprimepower, isperfectpower

# Levenshtein
export levenshtein, levenshtein!

include("abstract_types.jl")
include("rref.jl")
include("messages.jl") # implicitly exports distance.jl and powers.jl
include("bounds.jl")
include("algebra.jl")
include("levenshtein.jl")

end # end module

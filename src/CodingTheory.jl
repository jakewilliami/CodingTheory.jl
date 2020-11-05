#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
module CodingTheory

using Primes: isprime, primes
using LinearAlgebra: I
using Polynomials

include("utils.jl")

# Abstract types
export FinitePolynomial, AbstractCode, Alphabet, Messages, no_round, CodeUniverse, UniverseParameters

# RREF
export rref, rref!

# Messages, Distance, and Primes
export rate, sphere_covering_bound, sphere_packing_bound,
        construct_ham_matrix, isperfect, ishammingperfect, isgolayperfect,
        get_codewords_greedy, get_codewords_random, get_all_words,
        get_codewords
export hamming_distance, hamming_ball, code_distance, t_error_detecting,
        t_error_correcting, find_error_detection_max, find_error_correction_max
export isprimepower

# Algebra
export FinitePolynomial, list_polys, multiplication_table, list_span, islinear, isirreducible, normal_form!, normal_form, equivalent_code!, equivalent_code, generator!, generator, parity_check, syndrome, isincode

# Levenshtein
export levenshtein, levenshtein!

include("abstract_types.jl")
include("rref.jl")
include("messages.jl") # implicitly exports distance.jl and primes.jl
include("algebra.jl")
include("levenshtein.jl")

end # end module

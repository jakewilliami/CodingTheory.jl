#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
module CodingTheory

export FinitePolynomial, AbstractCode, Alphabet, Messages, no_round
export Alphabet, Messages, rate, sphere_covering_bound, sphere_packing_bound,
        construct_ham_matrix, isperfect, ishammingperfect, isgolayperfect,
        get_codewords_greedy, get_codewords_random, get_all_words,
        get_codewords
export FinitePolynomial, list_polys, multiplication_table, list_span, islinear,
        isirreducible
export hamming_distance, hamming_ball, code_distance, t_error_detecting,
        t_error_correcting, find_error_detection_max, find_error_correction_max
export levenshtein, levenshtein!
export normal_form!, normal_form, equivalent_code!, equivalent_code, generator!,
        generator, parity_check, syndrome, isincode
export rref, rref!

include(joinpath(dirname(@__FILE__), "abstract_types.jl"))
include(joinpath(dirname(@__FILE__), "messages.jl"))
include(joinpath(dirname(@__FILE__), "distance.jl"))
include(joinpath(dirname(@__FILE__), "fields.jl"))
include(joinpath(dirname(@__FILE__), "algebra.jl"))
include(joinpath(dirname(@__FILE__), "rref.jl"))

end # end module

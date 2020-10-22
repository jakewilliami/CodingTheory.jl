#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
module CodingTheory

export Alphabet, Messages, rate, sphere_covering_bound, sphere_packing_bound,
        construct_ham_matrix, isperfect, ishammingperfect, isgolayperfect,
        get_codewords, get_all_words
export list_polys, multiplication_table, list_span, islinear, isirreducible
export hamming_distance, hamming_ball, code_distance, t_error_detecting,
        t_error_correcting, find_error_detection_max, find_error_correction_max
export levenshtein, levenshtein!
export rref, rref!

include(joinpath(dirname(@__FILE__), "messages.jl"))
include(joinpath(dirname(@__FILE__), "distance.jl"))
include(joinpath(dirname(@__FILE__), "fields.jl"))
include(joinpath(dirname(@__FILE__), "rref.jl"))

end # end module

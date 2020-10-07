#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
module CodingTheory

export Alphabet, Messages, rate
export list_polys, multiplication_table, list_span, islinear, isirreducible
export hamming_distance, hamming_ball, code_distance, t_error_detecting, t_error_correcting
export levenshtein, levenshtein!

include(joinpath(dirname(@__FILE__), "messages.jl"))
include(joinpath(dirname(@__FILE__), "distance.jl"))
include(joinpath(dirname(@__FILE__), "fields.jl"))

end # end module

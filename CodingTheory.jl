#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
module CodingTheory

include(joinpath(dirname(@__FILE__), "distance.jl"))
include(joinpath(dirname(@__FILE__), "fields.jl"))

end # end module

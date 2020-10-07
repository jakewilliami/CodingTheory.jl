#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

@inline function _allequal_length_(A::AbstractArray)::Bool
    length(A) < 2 && return true
	
    e1 = A[1]
	
    @inbounds for i in 2:length(A)
        isequal(length(A[i]), length(e1)) || return false
    end
	
    return true
end

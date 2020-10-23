#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "utils.jl"))

#=
Gauss-Jordan elimination over finite fields, modulo n
=#

function rref!(A::Matrix{Int},
    n::Integer;
    colswap::Bool=false,
    verbose::Bool=false)::Matrix{Int}
    
    nrows, ncols = size(A)
    i = j = 1
    
    while i ≤ nrows && j ≤ ncols
        # Ignore zero rows.
        if isnothing(__findfirstnonzero(A[i, :]))
            i += 1
            continue
        end
        
        # Rule 1: Swap with the row above if out of order.
        if i > 1
            if isnothing(__findfirstnonzero(A[i - 1, :]))
                A[i,:], A[i - 1,:] = A[i - 1,:], A[i,:]
                if verbose
                    println("r$(i) ⟷ r$(i - 1)")
                end
                continue
            end
        end
        
        # Rule 2: Normalize each row
        s = __findfirstnonzero(A[i,:])
        α = invmod(A[i, s], n)
        A[i,:] .= mod.(A[i,:] * α, n)
        if verbose
            if ! iszero(α)
                isone(α) || println("r$(i) × $(α)")
            end
        end
        
        # Rule 3: Subtract it from the others
        s = __findfirstnonzero(A[i,:])
        for k in 1:nrows
            if i ≠ k
                β = A[k, s]
                A[k,:] .= mod.(A[k,:] - β * A[i,:], n)
                if verbose
                    if ! iszero(β)
                        if isone(β)
                            println("r$(k) - r$(i)")
                        else
                            println("r$(k) - $(β)r$(i)")
                        end
                    end
                end
            end
        end
        i += 1
        
        # Rule 4: Swap columns if needed (and allowed)
        if colswap
            if iszero(A[i, j])
                s = __findfirstnonzero(A[i,:])
                A[:,s], A[:,j] = A[:,j], A[:,s]
                if verbose
                    println("c$(j) ⟷ c$(s)")
                end
            end
        end
        
        # increment counters
        i += 1
        j += 1
    end
    
    if verbose
        println()
    end

    return A
end

rref(A::Matrix{Int},
    n::Integer;
    colswap::Bool=false,
    verbose::Bool=false
    )::Matrix{Int} = rref!(copy(A), n::Integer; colswap=colswap, verbose=verbose)

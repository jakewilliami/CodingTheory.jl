#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
using LinearAlgebra

findfirstnonzero(row::Vector)::Integer = findfirst(x -> ! iszero(x), row::Vector)

function rref!(A::Matrix{Int},
    n::Integer;
    colswap::Bool=false,
    verbose::Bool=false)::Matrix{Int}
    
    nrows, ncols = size(A)
    i = j = 1
    
    while i ≤ nrows && j ≤ ncols
        # Ignore zero rows.
        if isnothing(findfirstnonzero(A[i, :])) i += 1 end
        
        # Rule 1: Swap with the row above if out of order.
        if i > 1
            if isnothing(findfirstnonzero(A[i - 1, :]))
                A[i,:], A[j,:] = A[j,:], A[i,:]
                if verbose
                    println("r$(i) ⟷ r$(j)")
                end
            end
        end
        
        # Rule 2: Normalize each row
        s = findfirstnonzero(A[i,:])
        α = invmod(A[i, s], n)
        A[i,:] .= mod.(A[i,:] * α, n)
        if verbose
            isone(α) || println("r$(i) × $(α)")
        end
        
        # Rule 3: Subtract it from the others
        s = findfirstnonzero(A[i,:])
        for k in 1:nrows
            if i ≠ k
                β = A[k, s]
                A[k,:] .= mod.(A[k,:] - β * A[i,:], n)
                if verbose
                    if isone(β)
                        println("r$(k) - r$(i)")
                    else
                        println("r$(k) - $(β)r$(i)")
                    end
                end
            end
        end
        
        # Rule 4: Swap columns if needed (and allowed)
        if colswap
            if iszero(A[i, j])
                s = findfirstnonzero(A[i,:])
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

    return A
end

rref(A::Matrix{Int},
    n::Integer;
    colswap::Bool=false,
    verbose::Bool=false) = rref!(copy(A), n::Integer; colswap=colswap, verbose=verbose)

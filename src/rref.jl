#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "utils.jl"))

#=
Gauss-Jordan elimination over finite fields, modulo n
=#

function swaprows!(A::Matrix, i::Integer, j::Integer)
    for k in 1:size(A, 2)
        A[i, k], A[j, k] = A[j, k], A[i, k]
    end
    
    return nothing
end

function swapcols!(A::Matrix, i::Integer, j::Integer)
    for k in 1:size(A, 1)
        A[k, i], A[k, j] = A[k, j], A[k, i]
    end
    
    return nothing
end

function rref!(A::Matrix{Int},
    n::Integer;
    colswap::Bool=false,
    verbose::Bool=false,
    vverbose::Bool=false)::Matrix{Int}
    
    nrows, ncols = size(A)
    i = j = 1
    
    while i ≤ nrows && j ≤ ncols
        # Rule 1: Swap zero rows if out of order and ensure leading ones cascade down diagonally.
        s = findfirst(!iszero, A[i:end, :])
        isnothing(s) && break
        swaprows!(A, i, i + s[1] - 1)
        if verbose
            println("r$(i) ⟷  r$(i + s[1] - 1)")
            if vverbose
                displaymatrix(A); println()
            end
        end
        
        # Rule 2: Normalize each row
        s = findfirst(!iszero, A[i,:])
        isnothing(s) && break
        α = invmod(A[i, s], n)
        A[i,:] .= mod.(A[i,:] .* α, n)
        if verbose
            if ! iszero(α)
                isone(α) || println("r$(i) ⟵  r$(i) × $(α)")
                if vverbose && ! isone(α)
                    displaymatrix(A); println()
                end
            end
        end
        
        # Rule 3: Subtract it from the others
        s = findfirst(!iszero, A[i,:])
        isnothing(s) && break
        for k in 1:nrows
            if i ≠ k
                β = A[k, s]
                A[k, :] .= mod.(A[k, :] .- β .* A[i, :], n)
                if verbose
                    if ! iszero(β)
                        if isone(β)
                            println("r$(k) ⟵  r$(k) - r$(i)")
                        else
                            println("r$(k) ⟵  r$(k) - $(β)r$(i)")
                        end
                        if vverbose
                            displaymatrix(A); println()
                        end
                    end
                end
            end
        end
        
        # Rule 4: Swap columns if needed (and allowed)
        if colswap
            if iszero(A[i, j])
                s = findfirst(!iszero, A[i, :])
                swapcols!(A, j, s)
                if verbose
                    println("c$(j) ⟷  c$(s)")
                    if vverbose
                        displaymatrix(A); println()
                    end
                end
            end
        end
        
        # increment row and column counter respectively
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
    verbose::Bool=false,
    vverbose::Bool=false
    )::Matrix{Int} = rref!(copy(A), n::Integer; colswap=colswap, verbose=verbose, vverbose=vverbose)

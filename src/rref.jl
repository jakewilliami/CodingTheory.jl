#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

#=
Gauss-Jordan elimination over finite fields, modulo n
=#

function swaprows!(A::Matrix{T}, i::Int, j::Int) where T
    for k in axes(A, 2)
        A[i, k], A[j, k] = A[j, k], A[i, k]
    end
    
    return A
end

function swapcols!(A::Matrix{T}, i::Int, j::Int) where T
    for k in axes(A, 1)
        A[k, i], A[k, j] = A[k, j], A[k, i]
    end
    
    return A
end

"""
    rref!(A::Matrix{Int}, n::Int; colswap::Bool=false, verbose::Bool=false, vverbose::Bool=false) -> Matrix{Int}
    
Performs Gauss-Jordan Elimination on a matrix A.  *This directly changes the matrix A.  Use `rref` for a non-mutating version of this function.*

Parameters:
  - A::Matrix{Int}: A matrix of Ints you wish to perform Gauss-Jordan elimiation on.
  - n::Int: The modulus of the finite field you are working under.
  - colswap::Bool (kwarg): Whether or not you allow for column swapping.
  - verbose::Bool (kwarg): Print the row operations.
  - vverbose::Bool(kwarg): Print the intermediate matrices of the algorithm.
  
Returns:
  - Matrix{Int}: a matrix in row echelon form.
"""
function rref!(A::Matrix{Int},
    n::Int;
    colswap::Bool = false,
    verbose::Bool = false,
    vverbose::Bool = false)
    
    i = j = 1
    
    while i ≤ size(A, 1) && j ≤ size(A, 2)
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
        for k in axes(A, 1)
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

"""
    rref(A::Matrix{Int}, n::Int; colswap::Bool=false, verbose::Bool=false, vverbose::Bool=false) -> Matrix{Int}
    
Performs Gauss-Jordan Elimination on a matrix A.

Parameters:
  - A::Matrix{Int}: A matrix of Ints you wish to perform Gauss-Jordan elimiation on.
  - n::Int: The modulus of the finite field you are working under.
  - colswap::Bool (kwarg): Whether or not you allow for column swapping.
  - verbose::Bool (kwarg): Print the row operations.
  - vverbose::Bool(kwarg): Print the intermediate matrices of the algorithm.
  
Returns:
  - Matrix{Int}: a matrix in row echelon form.
"""
rref(A::Matrix{Int},
    n::Int;
    colswap::Bool = false,
    verbose::Bool = false,
    vverbose::Bool = false
    ) = rref!(copy(A), n::Int; colswap = colswap, verbose = verbose, vverbose = vverbose)

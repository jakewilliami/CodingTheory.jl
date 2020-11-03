#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

# Adapted from Levenshtein.jl

"""
    levenshtein(source::AbstractString, target::AbstractString) -> Integer
    levenshtein(source::AbstractString, target::AbstractString, cost::Real) -> Integer
    levenshtein(
        source::AbstractString,
        target::AbstractString,
        deletion_cost::R,
        insertion_cost::S,
        substitution_cost::T) -> Integer
    levenshtein!(
        source::AbstractString,
        target::AbstractString,
        deletion_cost::R,
        insertion_cost::S,
        substitution_cost::T,
        costs::Matrix = Array{promote_type(R, S, T)}(undef, 2, length(target) + 1)
    ) -> Integer
    
Computes the Levenshtein distance.  *These methods are adapted from Levenshtein.jl, by Roger Tu.*
"""
function levenshtein(source::AbstractString, target::AbstractString)
    return levenshtein(source, target, 1)
end

function levenshtein(source::AbstractString, target::AbstractString, cost::Real)
    return levenshtein(source, target, cost, cost, cost)
end

function levenshtein(
    source::AbstractString,
    target::AbstractString,
    deletion_cost::R,
    insertion_cost::S,
    substitution_cost::T) where {R <: Real, S <: Real, T <: Real}
    
    return levenshtein!(target, source, insertion_cost, deletion_cost, substitution_cost)
end

function levenshtein!(
    source::AbstractString,
    target::AbstractString,
    deletion_cost::R,
    insertion_cost::S,
    substitution_cost::T,
    costs::Matrix = Array{promote_type(R, S, T)}(undef, 2, length(target) + 1)
) where {R <: Real, S <: Real, T <: Real}

    cost_type = promote_type(R,S,T)
    
    if length(source) < length(target)
        # Space complexity of function = O(length(target))
        return levenshtein!(target, source, insertion_cost, deletion_cost, substitution_cost, costs)
    else
        if iszero(length(target))
            return length(source) * deletion_cost
        else
            old_cost_index = 1
            new_cost_index = 2

            costs[old_cost_index, 1] = 0
            for i in 1:length(target)
                costs[old_cost_index, i + 1] = i * insertion_cost
            end

            i = 0
            for r in source
                i += 1

                # Delete i characters from source to get empty target
                costs[new_cost_index, 1] = i * deletion_cost

                j = 0
                for c in target
                    j += 1

                    deletion = costs[old_cost_index, j + 1] + deletion_cost
                    insertion = costs[new_cost_index, j] + insertion_cost
                    substitution = costs[old_cost_index, j]
                    if r ≠ c
                        substitution += substitution_cost
                    end

                    costs[new_cost_index, j + 1] = min(deletion, insertion, substitution)
                end

                old_cost_index, new_cost_index = new_cost_index, old_cost_index
            end

            new_cost_index = old_cost_index
            
            return costs[new_cost_index, length(target) + 1]
        end
    end
end

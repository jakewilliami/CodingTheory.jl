#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
include(joinpath(dirname(@__FILE__), "distance.jl"))
include(joinpath(dirname(@__FILE__), "utils.jl"))

struct Alphabet
    Σ::Union{AbstractArray, AbstractString}
    
    function Alphabet(Σ::AbstractArray)
        try
            Σ = parse.(Int, Σ)
        catch
        end
        
        new(Σ)
    end # end constructor function
    
    function Alphabet(Σ::AbstractString)
        Σ = collect(Σ)
        
        try
            Σ = parse.(Int, Σ)
        catch
            Σ = string.(Σ)
        end
        
        new(Σ)
    end # end constructor function
end # end struct

struct Messages
    ℳ::AbstractArray
    length::Integer
    
    function Messages(ℳ::AbstractArray)
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(ℳ) || throw(error("$(message_length_error)"))
        
        length = length(ℳ[rand(1:length(ℳ))]) # choose arbitrary message in the list of messages
        
        new(ℳ, length)
    end # end constructor function
end

#=
n is the word length
q is the number of symbols in the code
M is the size/number of elements in the code
=#
function rate(q::Integer, M::Integer, n::Integer)
    return log(q, M) / n
end

function spheres(q::Integer, n::Integer, r::Integer)
    return sum([((q - 1)^i) * binomial(n, i) for i in 0:r])
end

function sphere_covering_bound(q::Integer, n::Integer, d::Integer)
    return Int(ceil((q^n) / spheres(q, n, d - 1)))
end

function sphere_packing_bound(q::Integer, n::Integer, d::Integer)
    return Int(floor((q^n) / spheres(q, n, Int(floor((d - 1) / 2)))))
end

function construct_ham_matrix(r::Integer, q::Integer)::Matrix
    ncols = Int(floor((q^r - 1) / (q - 1)))
    M = Matrix{Int}(undef, r, ncols)
    
    for i in 1:ncols
        M[:,i] = reverse(digits(parse(Int, string(i, base = q)), pad = r), dims = 1)
    end
    
    return M
end

#=
n is the word length
k is the dimension of the code (i.e., )
d is the distance of the code
If C ⊆ 𝔽_q^n, then C is over the field modulo q
=#
function isperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)::Bool
    
    M = q^k
    
    if isequal(sphere_packing_bound(q, n, d), M)
        return true
    else
        return false
    end
end

function ishammingperfect(r::Integer, q::Integer)::Bool
    
    n = 2^r - 1
    k = n - r
    M = q^k
    d = size(construct_ham_matrix(r, q))[1] # the number of rows of the hamming matrix (which is, by design, linearly independent)
    d = r
    
    println(n)
    println((q^r - 1) / (q - 1))
    
    if isequal(n, (q^r - 1) / (q - 1)) && isequal(d, 3) && isequal(M, q^(((q^r - 1) / (q - 1)) - r))
        return true
    end
    
    return false
end

function ishammingperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)::Bool
    
    M = q^k
    r = log(ℯ, ((n * log(ℯ, 1)) / (log(ℯ, 2))) + 1) / log(ℯ, 2)
    
    if isequal(n, (q^(r - 1)) / (q - 1)) && isequal(d, 3) && isequal(M, q^(((q^r - 1) / (q - 1)) - r))
        return true
    end
    
    return false
end

function isgolayperfect(
    n::Integer,
    k::Integer,
    d::Integer,
    q::Integer)::Bool
    
    M = q^k
    
    if isequal(q, 2) && isequal(n, 23) && isequal(d, 7) && isequal(M, 2^12)
        return true
    end
    
    if isequal(q, 3) && isequal(n, 11) && isequal(d, 5) && isequal(M, 3^6)
        return true
    end
    
    return false
end

function our_search()
    for n in 2:30, k in 2:30, d in 3:30
        if isperfect(n, k, d, 6)
            return n, k, d, 6
        end
    end
    
    return nothing
end

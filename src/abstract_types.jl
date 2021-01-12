#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

"""
```julia
abstract type FiniteField end
```
"""
abstract type FiniteField end

"""
```julia
abstract type AbstractCode end
```
"""
abstract type AbstractCode end

"""
```julia
NonStaticAbstractWord{N, T} = Union{NTuple{N, T}, AbstractVector{T}, AbstractString} where {N, T}
```
"""
NonStaticAbstractWord{N, T} = Union{NTuple{N, T}, AbstractVector{T}, AbstractString} where {N, T}

"""
``julia
mutable struct Word{N, T}
Word(w::NTuple{N, T})
Word(w::AbstractVector{T})
Word(w::AbstractString)
Word(i::T...)
```

A Word is an `StaticArrays.MVector` which is efficient (like tuple) yet mutable.
"""
mutable struct Word{N, T}
    w::MVector{N, T}
    
    Word(w::NTuple{N, T}) where {N, T} = new{N, T}(MVector{N, T}(w))
    Word(w::AbstractVector{T}) where {T} =
        (len = length(w); new{len, T}(MVector{len, T}(w)))
    Word(w::AbstractString) =
        MVector{length(w), eltype(w)}(collect(w))
    Word(i::T...) where {T} =
        (len = length(i); new{len, T}(MVector{len, T}(i)))
end

# Indexing Interface
Base.getindex(W::Word{N, T}, i::Int) where {N, T} = getindex(W.w, i)
Base.setindex!(W::Word{N, T}, v, i::Int) where {N, T} = setindex!(W.w, v, i)
Base.firstindex(W::Word{N, T}) where {N, T} = firstindex(W.w)
Base.lastindex(W::Word{N, T}) where {N, T} = lastindex(W.w)

# Abstract Array Interface
Base.size(::Word{N, T}) where {N, T} = tuple(N)
Base.length(::Word{N, T}) where {N, T} = N

"""
```julia
AbstractWord{N} = Union{NonStaticAbstractWord{N, T}, MVector{T}} where {T}
```
"""
AbstractWord{N, T} = Union{Word{N, T}, NonStaticAbstractWord{N, T}, MVector{N, T}} where {N, T}

isword(w) = w isa Word
isword(i...) = isword(i)
isabstractword(w) = w isa AbstractWord
isabstractword(i...) = isabstractword(i)

"""
```julia
Codewords{N} <: AbstractCode
```

Simply a wrapper type for a vector of abstract words of length `N`.
"""
Codewords{N} = Vector{AbstractWord{N, T}} where T
# Codewords{N} = Vector{Word{N, T}} where T


"""
```julia
struct Alphabet <: AbstractCode
```

Has the parameter Σ, which is the alphabet; a collection of strings, characters, symbols, or Ints.

---

```julia
Alphabet(Σ::AbstractArray)
```

A constructor method for the struct Alphabet.  Takes an array of letters in the alphabet, and attempts to parse them as 64-bit Ints.
    
---

```julia
Alphabet{N}(Σ::AbstractArray)
Alphabet{N}(Σ::AbstractString)
```

A constructor method for the struct Alphabet.  Takes in a symbols and splits it into constituent characters.  Those symbols are the letters in the alphabet.  `N` is the number of letters in the alphabet
"""
struct Alphabet{N} <: AbstractVector{Symbol} where N
    Σ::AbstractVector{Symbol}

    Alphabet(Σ::Union{Vector{T}, String}) where {T} =
        (Σ_unique = Set(Σ); new{length(Σ_unique)}(ensure_symbolic(Σ_unique)))
end # end struct

# Indexing Interface
Base.getindex(A::Alphabet, i::Int) = getindex(A.Σ, i)
Base.setindex!(A::Alphabet, v, i::Int) = setindex!(A.Σ, v, i)
Base.firstindex(A::Alphabet) = firstindex(A.Σ)
Base.lastindex(A::Alphabet) = lastindex(A.Σ)

# Abstract Array Interface
Base.size(A::Alphabet{N}) where {N} = tuple(N)
Base.length(A::Alphabet{N}) where {N} = N
Base.rand(A::Alphabet{N}) where {N} = rand(A.Σ)

"""
```julia
gensym(q::Int) -> Vector{Symbol}
```

Generates a vector of unique symbols of length q.
"""
Base.gensym(q::Int) = Symbol[gensym() for _ in 1:q]

"""
```julia
genalphabet(q::Int)
```

Generates an alphabet of q unique symbols.
"""
genalphabet(q::Int) = Alphabet(gensym(q))

"""
```julia
struct UniverseParameters <: AbstractCode
```

Defines a structure for the messages in the code.  Parameters are the alphabet `Σ`, size of alphabet `q`, and block length `n`

```julia
UniverseParameters(Σ::Alphabet, n::Int)
UniverseParameters(Σ::AbstractArray, n::Int)
UniverseParameters(Σ::Alphabet, q::Int, n::Int)
UniverseParameters(Σ::AbstractArray, q::Int, n::Int)
UniverseParameters(q::Int, n::Int)
```

An inner constructor function on the structure `UniverseParameters`.
"""
struct UniverseParameters <: AbstractCode
    Σ::Alphabet{N} where N
    q::Int # size of alphabet
    n::Int # block length
    
    UniverseParameters(Σ::Alphabet{N}, n::Int) where {N} = new(Σ, N, n)
    UniverseParameters(Σ::Alphabet{N}, q::Int, n::Int) where {N} = new(Σ, q, n)
    UniverseParameters(Σ::AbstractVector{T}, n::Int) where {T} = (Σ = Alphabet(Σ); new(Σ, length(Σ), n))
    UniverseParameters(Σ::AbstractVector{T}, q::Int, n::Int) where {T} = new(Σ, q, n)
    UniverseParameters(q::Int, n::Int) = new(genalphabet(q), q, n)
end

# Iteration inteface functions
Base.length(iter::UniverseParameters) = big(iter.q)^iter.n
Base.eltype(iter::UniverseParameters) = NTuple{iter.n, Symbol}

# Other Base interface functions for UniverseParameters
Base.rand(𝒰::UniverseParameters) = ntuple(_ -> rand(𝒰.Σ), 𝒰.n)

"""
```julia
rand(𝒰::UniverseParameters, C::AbstractArray) -> Tuple
```

Given universe parameters 𝒰 and a code C, return a tuple including
  - A random word in C;
  - A random letter in the alphabet; and
  - A random index in the block length.
"""
Base.rand(𝒰::UniverseParameters, C::AbstractArray) = tuple(Word(rand(C)), rand(𝒰.Σ), rand(1:𝒰.n))

"""
```julia
struct CodeUniverseIterator <: AbstractCode
```

A structure used to iterate through a code universe, with specified universe parameters.  E.g.,

```julia
for c in CodeUniverseIterator(["a", "b", "c"], 3)
    println(c)
end
```

Fields:
  - `𝒰::UniverseParameters`

Methods:

```julia
CodeUniverseIterator(𝒰::UniverseParameters)
CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, q::Int, n::Int)
CodeUniverseIterator(Σ::Union{Alphabet, AbstractArray}, n::Int)
CodeUniverseIterator(q::Int, n::Int)
```
"""
struct CodeUniverseIterator <: AbstractCode
    𝒰::UniverseParameters

    function CodeUniverseIterator(𝒰::UniverseParameters)
        Σ = 𝒰.Σ
        return Iterators.product(Vector{eltype(Σ)}[Σ for _ in 1:𝒰.n]...)
    end

    CodeUniverseIterator(Σ::Union{Alphabet{N}, AbstractVector{T}}, q::Int, n::Int) where {N, T} =
        CodeUniverseIterator(UniverseParameters(Σ, q, n))
    CodeUniverseIterator(Σ::Union{Alphabet{N}, AbstractVector{T}}, n::Int) where {N, T} =
        CodeUniverseIterator(UniverseParameters(Σ, n))
    CodeUniverseIterator(q::Int, n::Int) =
        CodeUniverseIterator(UniverseParameters(q, n))
end

"""
```julia
struct CodeUniverse <: AbstractCode
```
    
Defines a structure for the messages in the code.  Parameters are the abstract array of messages `𝒰`, and the length of the messages `n`.

```julia
CodeUniverse(𝒰::AbstractArray, Σ::Alphabet)
```
    
An inner constructor function on the structure `CodeUniverse`.
"""
struct CodeUniverse <: AbstractCode
    𝒰::Codewords{M} where {M}
    Σ::Alphabet{N} where N
    q::Int # alphabet size
    n::Int # block length
    
    function CodeUniverse(𝒰::Codewords{M}, Σ::Alphabet{N}) where {N, M}
        message_length_error = "We have fixed the block length of each message.  Please ensure all messages are of equal length."
        _allequal_length_(𝒰) || throw(error(message_length_error))
    
        new(𝒰, Σ, length(Σ), length(first(𝒰)))
    end # end constructor function
end

"""
```julia
struct Rounding end
```

An abstract structure used as a parameter in a non-rounding method of `hamming_bound`.  Use `no_round` for this; a predefiend variable of the structure `Rounding`.
"""
struct Rounding end
no_round = Rounding()

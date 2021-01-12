#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

"""
```julia
displaymatrix(M::AbstractArray)
```
	
Displays a matrix `M` in a compact form from the terminal.
"""
function displaymatrix(M::AbstractArray)
    return show(IOContext(stdout, :limit => true, :compact => true, :short => true), "text/plain", M); print("\n")
end

_Iterable{T} = Union{AbstractArray{T}, NTuple{N, T}} where N

"""
```julia
allequal_length(A) -> Bool
allequal_length(a, b...) -> Bool
```

Check that all elements in a list are of equal length.
"""
@inline function allequal_length(A::_Iterable{T}) where T
    length(A) < 2 && return true
	
    @inbounds for i in 2:length(A)
        isequal(length(first(A)), length(A[i])) || return false
    end
	
    return true
end
@inline allequal_length(a::T...) where {T} = allequal_length(a)

"""
```julia
allequal(A) -> Bool
allequal(a, b...) -> Bool
```

Check that all elements in a list are equal to each other.
"""
@inline function allequal(A::_Iterable{T}) where T
    length(A) < 2 && return true
    
    @inbounds for i in 2:length(A)
        first(A) ≠ A[i] && return false
    end
    
    return true
end
@inline allequal(a::T...) where {T} = allequal(a)

"""
```julia
aredistinct(A) -> Bool
aredistinct(a, b...) -> Bool
```

Check that all elements in a list are distinct from every other element in the list.
"""
@inline function aredistinct(A::_Iterable{T}) where T
    length(A) < 2 && return true
    
    while ! iszero(length(A))
        a = pop!(A)
        a ∉ A || return false
    end
    
    return true
end
@inline aredistinct(a::T...) where {T} = aredistinct(a)

"""
```julia
arelessthan(x::Number, A) -> Bool
arelessthan(x::Number, a, b...) -> Bool
```

Check that all elements in a list are less than a given x.
"""
@inline function arelessthan(x::Number, A::_Iterable{T}) where T
    @inbounds for a in A
        a < x || return false
    end
    
    return true
end
@inline arelessthan(x::Number, a::Number...) = arelessthan(x, a)

"""
```julia
areequalto(x::Number, A::AbstractArray) -> Bool
areequalto(x::Number, a, b...) -> Bool
```

Check that all elements in a list are equal to a given x.
"""
@inline function areequalto(x, A::_Iterable{T}) where T
    @inbounds for a in A
        a != x || return false
    end
    
    return true
end
@inline areequalto(x, a::Number...) = areequalto(x, a)

"""
```julia
deepsym(a::AbstractArray)
```

Convert inner-most elements into symbols
"""
deepsym(a) = Symbol.(a)
deepsym(a::_Iterable{T}) where {T} = deepsym.(a)

"""
```julia
deepeltype(a::AbstractArray) -> Type
```

Returns the type of the inner-most element in a nested array structure.
"""
deepeltype(a) = deepeltype(typeof(a))
deepeltype(::Type{T}) where {T <: _Iterable{R}} where {R} = deepeltype(eltype(T))
deepeltype(::Type{T}) where T = T

"""
```julia
ensure_symbolic!(Σ) -> typeof(Σ)
```

Ensures that the inner-most elements of a nested array structure are of the type `Symbol`.  *This is a mutating function.  Use its twin, non-mutating function, `ensure_symbolic`, if you need a non-mutating version of this.*
"""
ensure_symbolic!(Σ) = deepeltype(Σ) isa Symbol ? Σ : deepsym(Σ)

"""
```julia
ensure_symbolic(Σ) -> typeof(Σ)
```

Ensures that the inner-most elements of a nested array structure are of the type `Symbol`.
"""
ensure_symbolic(Σ) = ensure_symbolic!(copy(Σ))

"""
```julia
_has_identity(M::Matrix, n::Union{T, Base.OneTo{T}, AbstractRange{T}}) where {T <: Integer}
```
Inner function on `has_identity` function.  Will return a tuple of:
  - Whether or not the matrix has an identity;
  - Where that identity matrix starts; and
  - How big the identity matrix is.
"""
function _has_identity(M::Matrix, n::Union{T, Base.OneTo{T}, AbstractRange{T}}) where {T <: Integer}
	for ID_size in n # iterate through possible identity sizes; try identity matrixes from the maximum size down
		n isa Integer || (isone(ID_size) && continue) # skip I(1) unless that is the integer specifically given
		# println("trying identity of size $ID_size")
		max_starting_row, max_starting_col = size(M) .- (ID_size - 1)
		for i in CartesianIndices(M) # iterate through all possible starting positions
			# println("Checking position $i; is it in the bounds of ($max_starting_row, $max_starting_col)?")
			_row, _col = Tuple(i)
			if _row ≤ max_starting_row && _col ≤ max_starting_col # ensure the matrix at position i has enough space for an identity
				# println("Yes it is; now we need to check if $(M[_row:(_row + ID_size - 1), _col:(_col + ID_size - 1)]) is an identity matrix")
				isequal(M[_row:(_row + ID_size - 1), _col:(_col + ID_size - 1)], I(ID_size)) && return true, (_row, _col), ID_size
			end
		end
	end
	
	return false, ntuple(_ -> nothing, ndims(M)), nothing
end

"""
```julia
has_identity(M::Matrix) -> Bool
has_identity(M::Matrix, n::Integer) -> Bool
```

Checks if a matrix has an identity in it.  If given a number `n`, the function will specifically check if it has an identity matrix _of size n_ in `M`.
	
See `CodingTheory._has_identity` for more information.

---

### Examples

```julia
julia> A = [1 0 0 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]
3×6 Array{Int64,2}:
 1  0  0  2  3  0
 0  1  0  1  2  2
 0  0  1  4  3  0

julia> B = [1 0 1 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]
3×6 Array{Int64,2}:
 1  0  1  2  3  0
 0  1  0  1  2  2
 0  0  1  4  3  0

julia> C = [-96 -66 20 1 0 0; -65 59 -82 0 1 0; -16 87 -113 0 0 1]
3×6 Array{Int64,2}:
-96  -66    20  1  0  0
-65   59   -82  0  1  0
-16   87  -113  0  0  1

julia> D = [78 -99 125 -123 -111 -71 17; -115 78 40 -88 81 -40 78; -99 126 -54 1 0 0 24; -55 88 42 0 1 0 -8; 119 55 2 0 0 1 -92; -40 -21 -89 -79 59 -44 9]
6×7 Array{Int64,2}:
   78  -99  125  -123  -111  -71   17
 -115   78   40   -88    81  -40   78
  -99  126  -54     1     0    0   24
  -55   88   42     0     1    0   -8
  119   55    2     0     0    1  -92
  -40  -21  -89   -79    59  -44    9

julia> has_identity(A)
true

julia> has_identity(B)
true

julia> has_identity(B, 3) # no identity matrix of size 3 exists
false

julia> has_identity(C)
true

julia> has_identity(D)
true

julia> has_identity(D, 4)
false
```
"""
has_identity(M::Matrix) =
	_has_identity(M, reverse(minimum(axes(M))))[1]
has_identity(M::Matrix, n::T) where {T <: Integer} =
	_has_identity(M, n)[1]

"""
```julia
has_identity_on_left(M::Matrix) -> Bool
```

Checks if the left-hand side of a matrix contains an identity matrix.

---

### Examples

```julia
julia> A = [1 0 0 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]
3×6 Array{Int64,2}:
 1  0  0  2  3  0
 0  1  0  1  2  2
 0  0  1  4  3  0

julia> B = [-96 -66 20 1 0 0; -65 59 -82 0 1 0; -16 87 -113 0 0 1]
3×6 Array{Int64,2}:
 -96  -66    20  1  0  0
 -65   59   -82  0  1  0
 -16   87  -113  0  0  1

julia> has_identity_on_left(A)
true

julia> has_identity_on_left(B)
false
```
"""
has_identity_on_left(M::Matrix) = isequal(M[:, axes(M, 1)], I(size(M, 1))) ? true : false

"""
```julia
sizeof_perfect_code(q::Number, n::Number, d::Number) -> Number
```

Calculates the number of gigabytes required to store a perfect code of parameters q, n, and d.
"""
function sizeof_perfect_code(q::Int, n::Int, d::Int)
	return (sizeof(ntuple(_ -> gensym(), n)) * hamming_bound(big.([q, n, d])...)) / (2^30)
end
sizeof_perfect_code(q::Number, n::Number, d::Number) = sizeof_perfect_code(round.(BigInt, [q, n, d])...)

"""
```julia
sizeof_all_words(q::Number, n::Number) -> Number
```

Calculates the number of gigabytes required to store all unique words of length n from an alphabet of size q.
"""
function sizeof_all_words(q::Int, n::Int)
	return (sizeof(ntuple(_ -> gensym(), n)) * big(q)^n) / (2^30)
end
sizeof_all_words(q::Number, n::Number) = sizeof_all_words(round.(BigInt, (q, n))...)

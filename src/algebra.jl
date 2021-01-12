#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

"""
```julia
struct FinitePolynomial <: FiniteField
```

Has parameters `p`, which is an abstract polynomial, and `n` which is the modulus of the field under which the molynomial is defined.

---

```julia
FinitePolynomial(p::AbstractPolynomial, n::Int)
```

A constructor method for `FinitePolynomial`.  Takes in a polynomial `p` and a number `n`, and constructs a polynomial under modulo n.
"""
struct FinitePolynomial <: FiniteField
	p::AbstractPolynomial
	n::Int
	
	function FinitePolynomial(p::AbstractPolynomial, n::Int)
		p = Polynomial(mod.(p.coeffs, n))
		new(p, n)
	end
end

"""
```julia
mod(p::Polynomial, n::Int) -> Polynomial
```
	
Uses the `FinitePolynomial` constructor to return a polynomial `p` under modulus `n`.

Parameters:
  - `p::Polynomial`: The input polynomial.
  - `n::Int`: The modulus of the field.
  
Returns
  - Polynomial: A polynomial modulo n.
"""
Base.mod(p::Polynomial, n::Int) = FinitePolynomial(p, n).p

"""
```julia
Polynomial(A::Union{NTuple{N, T}, Vector{T}}, n::Int) -> Polynomial
```
	
Constructs a polynomial under modulo `n`.

Parameters:
  - `A::Union{Tuple, AbstractArray}`: The polynomial coefficients.
  - `n::Int`: The modulus of the field.
  
Returns
  - `Polynomial`: A polynomial modulo n.
"""
Polynomial(A::Union{NTuple{N, T}, Vector{T}}, n::Int) where {N, T} = mod(Polynomial(A), n)

"""
```julia
list_polys(n::Int, m::Int) -> Array
```
	
Lists all polynomials of degree less than to `n` under modulo `m`.

Parameters:
  - `n::Int`: Highest degree of polynomial.
  - `m::Int`: The modulus of the field.
  
Returns:
  - `Array`: An array of polynomials of degree less than n, under modulo m.
"""
function list_polys(n::Int, m::Int)
	return collect(Polynomial(collect(t)) for t in Iterators.product([0:(m-1) for i in 1:n]...))
end

"""
```julia
multiplication_table(degree::Int, modulo::Int) -> Matrix
```
	
Returns a table (matrix) of the multiplication of all combinations of polynomials for degree less than `degree`, under modulo `modulo`.

Parameters:
  - `degree::Int`: Highest degree of polynomial.
  - `modulo::Int`: The modulus of the field.
  
Returns:
  - `Matrix`: A multiplication table of all polynomials with degree less than n, under modulus.

---

### Examples

```julia
julia> julia> multiplication_table(2, 3) # multiplication table of all polynomials of degree less than 3 modulo 2
9×9 Array{Polynomial,2}:
 Polynomial(0)  Polynomial(0)        Polynomial(0)        …  Polynomial(0)                Polynomial(0)
 Polynomial(0)  Polynomial(1)        Polynomial(2)           Polynomial(1 + 2*x)          Polynomial(2 + 2*x)
 Polynomial(0)  Polynomial(2)        Polynomial(1)           Polynomial(2 + x)            Polynomial(1 + x)
 Polynomial(0)  Polynomial(x)        Polynomial(2*x)         Polynomial(x + 2*x^2)        Polynomial(2*x + 2*x^2)
 Polynomial(0)  Polynomial(1 + x)    Polynomial(2 + 2*x)     Polynomial(1 + 2*x^2)        Polynomial(2 + x + 2*x^2)
 Polynomial(0)  Polynomial(2 + x)    Polynomial(1 + 2*x)  …  Polynomial(2 + 2*x + 2*x^2)  Polynomial(1 + 2*x^2)
 Polynomial(0)  Polynomial(2*x)      Polynomial(x)           Polynomial(2*x + x^2)        Polynomial(x + x^2)
 Polynomial(0)  Polynomial(1 + 2*x)  Polynomial(2 + x)       Polynomial(1 + x + x^2)      Polynomial(2 + x^2)
 Polynomial(0)  Polynomial(2 + 2*x)  Polynomial(1 + x)       Polynomial(2 + x^2)          Polynomial(1 + 2*x + x^2)
```
"""
function multiplication_table(degree::Int, modulo::Int)
	polys = list_polys(degree, modulo)
	number_of_polys = length(polys)
	poly_matrix = Matrix{Polynomial}(undef, number_of_polys, number_of_polys)
	
	for i in 1:number_of_polys, j in 1:number_of_polys
		poly_matrix[i,j] = mod(polys[i] * polys[j], modulo)
	end

	return poly_matrix
end

function __list_span_inner(modulo::Int, u̲::Vector{T}...) where T <: Number
    n_vec = length(u̲)
    span = Vector{T}[zero(u̲[1])]
	
    for n in 1:n_vec
        new_span = copy(span)
        for v in span, λ in 1:(modulo - 1)
            w̲ = mod.(v + λ * u̲[n], modulo)
            if w̲ ∉ new_span
                push!(new_span, w̲)
            end
        end
        span = new_span
    end
	
	return span
end

"""
```julia
list_span(u̲::Vector, v̲::Vector,... modulo::Int) -> Array
```

Given any number of vectors `u̲`, `v̲`,... , prints all linear combinations of those vectors, modulo `modulo`.

Parameters:
  - `u̲::Vector`: One vector.
  - `v̲::Vector`: Another vector.
  - `...`: Other vectors.
  - `modulo::Int`: The modulus of the field.
  
Returns:
  - `Array`: All vectors in the span of u̲ and v̲, under modulo.

---

### Examples

```julia
julia> list_span([2, 1, 1], [1, 1, 1], 3) # list the span of two vectors modulo 3
9-element Array{Array{T,1} where T,1}:
 [0, 0, 0]
 [1, 1, 1]
 [2, 2, 2]
 [2, 1, 1]
 [0, 2, 2]
 [1, 0, 0]
 [1, 2, 2]
 [2, 0, 0]
 [0, 1, 1]
```
"""
list_span(args...) = __list_span_inner(last(args), args[1:end-1]...)

"""
```julia
islinear(C::Vector, modulo::Int; verbose::Bool = false) -> Bool
```
	
Determines whether a code `C` is a linear code (i.e., if it is closed under addition, scalar multiplication, and has the zero vector in it).

Parameters:
  - `C::Vector`: A code, typically consisting of multiple vectors or strings.
  - `modulo::Int`: The modulus of the field under which you are working.
  - `verbose::Bool` (kwarg): print the point at which C fails, if it does.
  
Returns:
  - `Bool`: Whether or not the code `C` is linear (true or false).

---

### Examples

```julia

julia> islinear([[0,0,0],[1,1,1],[1,0,1],[1,1,0]], 2) # checks whether a vector of vectors is linear/a subspace (modulo 2)
false
```
"""
function islinear(C::Vector, modulo::Int; verbose::Bool = false)
	allequal_length(C) || return false # not all codes are of the same length
	block_length = length(C[1])
	𝟎 = fill(0, block_length)
		
	if 𝟎 ∉ C
		verbose && println("The zero vector 0̲ is not in C.\n")
		return false # the zero vector is not in the code
	end
	
	for c̲ ∈ C
		for λ in 0:modulo-1
			if mod.(λ*c̲, modulo) ∉ C
				verbose && println(λ, " ⋅ ", c̲, " = ", mod.(λ*c̲, modulo), " ∉ C\n")
				return false # this code isn't closed under scalar multiplication
			end
		end
		
		for c̲ ∈ C, c̲′ ∈ C
			if c̲ ≠ c̲′
				if mod.(c̲ + c̲′, modulo) ∉ C
					verbose && println(c̲, " + ", c̲′, " = ", mod.(c̲ + c̲′, modulo), " ∉ C\n")
					return false # this code isn't closed under addition
				end
			end
		end
	end
	
	verbose && println()
	
	return true
end

"""
```julia
isirreducible(f::AbstractPolynomial, modulo::Int) -> Bool
```
	
Checks if a polynomial is irreducible.
	
Parameters:
  - `f::Polynomial`: The polynomial you need to check.
  - `modulo::Int`: The modulus under which you are working.
  
Returns:
  - `Bool`: Whether or not the polynomial is irreducible (true or false).

---

### Examples

```julia
julia> isirreducible(Polynomial([1, 1, 0, 0, 1]), 2) # is 1 + x + x^4 mod 2 irreducible?
true
```
"""
function isirreducible(f::Polynomial, modulo::Int)
	deg = length(f.coeffs) - 1
	f = mod(f, deg)
	polys = list_polys(deg, modulo)

	for a in polys, b in polys
		isequal(f, mod(a*b, modulo)) && return false
	end
	
	return true
end

"""
```julia
normal_form!(M::AbstractArray{Int}, n::Int) -> Matrix{Int}
```

Convert a matrix M into normal form under modulo n via Gauss-Jordan elimination.  *This directly changes the matrix M.  Use `normal_form` for a non-mutating version of this function.*

Parameters:
  - `M::AbstractArray{Int}`: A matrix of Ints.
  - `n::Int`: The modulus of the finite field.
  
Returns:
  - `Matrix{Int}`: A matrix in normal form from Gauss-Jordan elimination.

---

### Examples

```julia

julia> normal_form!([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) # computes rref colswap = false
4×7 Array{Int64,2}:
1  0  0  0  2  2  2
0  1  0  0  2  0  1
0  0  1  0  1  0  2
0  0  0  1  2  2  1
```
"""
function normal_form!(M::AbstractArray{Int}, n::Int)
	return rref!(M, n, colswap = false)
end

"""
```julia
normal_form(M::AbstractArray{Int}, n::Int) -> Matrix{Int}
```

Convert a matrix M into normal form under modulo n via Gauss-Jordan elimination.

Parameters:
  - `M::AbstractArray{Int}`: A matrix of Ints.
  - `n::Int`: The modulus of the finite field.
  
Returns:
  - `Matrix{Int}`: A matrix in normal form from Gauss-Jordan elimination.

---

### Examples

```julia

julia> normal_form([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) # computes rref colswap = false
4×7 Array{Int64,2}:
1  0  0  0  2  2  2
0  1  0  0  2  0  1
0  0  1  0  1  0  2
0  0  0  1  2  2  1
```
"""
function normal_form(M::AbstractArray{Int}, n::Int)
	normal_form!(copy(M), n)
end

"""
```julia
equivalent_code!(M::AbstractArray{Int}, n::Int) -> Matrix{Int}
```

Peforms Gauss-Jordan elimination on a matrix M, but allows for column swapping.  This constructs an "equivalent" matrix.  *This directly changes the matrix M.  Use `equivalent_code` for a non-mutating version of this function.*

Parameters:
  - `M::AbstractArray{Int}`: A matrix of Ints.
  - `n::Int`: The modulus of the finite field.
  
Returns:
  - `Matrix{Int}`: A which represents an "equivalent" code to that of the matrix M.


---

### Examples

```julia
julia> equivalent_code!([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) # computes rref colswap = true
4×7 Array{Int64,2}:
1  0  0  0  2  2  2
0  1  0  0  2  0  1
0  0  1  0  1  0  2
0  0  0  1  2  2  1
```
"""
function equivalent_code!(M::AbstractArray{Int}, n::Int)
	return rref!(M, n, colswap=true)
end

"""
```julia
equivalent_code(M::AbstractArray{Int}, n::Int) -> Matrix{Int}
```

Peforms Gauss-Jordan elimination on a matrix M, but allows for column swapping.

Parameters:
  - `M::AbstractArray{Int}`: A matrix of Ints.
  - `n::Int`: The modulus of the finite field.
  
Returns:
  - `Matrix{Int}`: A which represents an "equivalent" code to that of the matrix M.

---

### Examples

```julia
julia> equivalent_code([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) # computes rref colswap = true
4×7 Array{Int64,2}:
1  0  0  0  2  2  2
0  1  0  0  2  0  1
0  0  1  0  1  0  2
0  0  0  1  2  2  1
```
"""
function equivalent_code(M::AbstractArray{Int}, n::Int)
	equivalent_code!(copy(M), n)
end

"""
```julia
generator!(M::AbstractArray{Int}, n::Int; colswap::Bool = true) -> Matrix{Int}
```

Constructs a generator matrix of the code, depending on if you allow for column swapping or not.   This function uses `normal_form!` or `equivalent_code!`.  *This directly changes the matrix M.  Use `generator` for a non-mutating version of this function.*

Parameters:
  - `M::AbstractArray{Int}`: A matrix of Ints.
  - `n::Int`: The modulus of the finite field.
  - `colswap::Bool` (kwarg): A boolean flag indicating whether or not you allow for swapping of columns when constructing the generating matrix.
	
Returns:
  - `Matrix{Int}`: A generating matrix.
"""
function generator!(M::AbstractArray{Int}, n::Int; colswap::Bool = false)
	return ifelse(colswap, equivalent_code!(M, n), normal_form!(M, n))
end

"""
```julia
generator(M::AbstractArray{Int}, n::Int; colswap::Bool = true) -> Matrix{Int}
```

Constructs a generator matrix of the code, depending on if you allow for column swapping or not.  This function uses `normal_form` or `equivalent_code`.

Parameters:
  -` M::AbstractArray{Int}`: A matrix of Ints.
  - `n::Int`: The modulus of the finite field.
  - `colswap::Bool` (kwarg): A boolean flag indicating whether or not you allow for swapping of columns when constructing the generating matrix.
	
Returns:
  - `Matrix{Int}`: A generating matrix.
"""
function generator(M::AbstractArray{Int}, n::Int; colswap::Bool = true)
	generator!(copy(M), n, colswap=colswap)
end

"""
```julia
parity_check(M::AbstractArray{Int}, n::Int) -> Matrix{Int}
```
	
Constructs a parity check matrix.  This is calculated from taking the non-identity part of a matrix in normal form (or equivalent &mdash; see `generator`), transposing it, multiplying it by negative one, and appending to it an appropriate sized identity matrix.
	
Parameters:
  - `M::AbstractArray{Int}`: A matrix of Ints.
  - `n::Int`: The modulus of the finite field.
	
Returns:
  - `Matrix{Int}`: A parity check matrix.
"""
function parity_check(M::AbstractArray{Int}, n::Int)
	has_identity_on_left(M) || throw(error("This matrix is not in normal form.  Use normal_form or equivalent_code."))
	
	minus_Dᵀ = mod.(-one(Int) .* transpose(M[:, (size(M, 1) + 1):end]), n)
	return H = Int[minus_Dᵀ I(size(minus_Dᵀ, 1))]
end

"""
```julia
syndrome(v̲::Vector, Hᵀ::AbstractArray{Int}, n::Int) -> Matrix{Int}
```
	
Calculates the syndrome of a given word v̲ and a parity check matrix, transposed (Hᵀ), under modulo n.

Parameters:
  - `v̲::Vector`: A word in the code.
  - `Hᵀ::AbstractArray{Int}`: The transpose of a parity check matrix.
  
Returns:
  - `Vector`: The syndrome of a word in the code.
 
---

### Examples

```julia
julia> syndrome([0, 2, 1, 2, 0, 1, 0], [1 1 1; 1 0 2; 2 0 1; 1 1 2; 1 0 0; 0 1 0; 0 0 1], 3)
1×3 Array{Int64,2}:
0  0  0

julia> syndrome([0, 2, 1, 2, 0, 1, 0], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3)
1×3 Array{Int64,2}:
 0  0  0
```
"""
function syndrome(v̲::Vector, Hᵀ::AbstractArray{Int}, n::Int)
	return mod.(v̲' * Hᵀ, n)
end

"""
```julia
isincode(v̲::Vector, Hᵀ::AbstractArray{Int}, n::Int) -> Bool
```

If the syndrome of a code is the zero vector, then the word used to calculate the syndrome is in the code.

Parameters:
  - `v̲::Vector`: A word.
  - `Hᵀ::AbstractArray{Int}`: The transpose of a parity check matrix.
  
Returns:
  - `Bool`: If the word is in the code or not (true or false).

---

### Examples

```julia
julia> isincode([0, 2, 1, 2, 0, 1, 0], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3) # tests if the syndrome is equal to the zero vector, and is thus in the code
true
```
"""
function isincode(v̲::Vector, Hᵀ::AbstractArray{Int}, n::Int)
	return iszero(syndrome(v̲, Hᵀ, n))
end

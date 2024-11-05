# Perfect power function

@doc raw"""
```julia
isperfectpower(n::Integer) -> Bool
```

Given an integer `n`, returns `true` or `false` depending on whether or not it is a perfect power.

Here, a perfect power is some number of the form ``a^b``, where ``a, b \in \mathbb{N}``, and ``b > 1``.

This function is a wrapper around [Hecke.jl's really excellent and efficient `ispower` function](https://github.com/thofma/Hecke.jl/blob/master/src/Misc/Integer.jl#L413-L443).

!!! note
	
	It should be noted that there is another defintion for "perfect powers":
		> some number of the form ``p^b``, where ``p, b \in \mathbb{N}``, *and ``p`` is a prime number``*.
	Here, we call numbers of this form *prime powers*, or *proper perfect powers*.  By this definition, 36 is *not* a perfect power, because 6 is not prime.

See also: `CodingTheory.isperfectpower`.
"""
function isperfectpower(n::Integer)
    (n == 0) && return false
    e = first(ispower(n))
    return (e != 1 && e != 0)
end

@doc raw"""
```julia
isprimepower(n::Integer) -> Bool
```

Given an integer `n`, returns `true` or `false` depending on whether or not it is a prime power.

A prime power is some number of the form ``p^b``, where ``p, b \in \mathbb{N}``, and ``p`` is a prime number``.

This function is a wrapper around [Hecke.jl's really excellent and effcient `ispower` function](https://github.com/thofma/Hecke.jl/blob/master/src/Misc/Integer.jl#L756-L769).

See also: `CodingTheory.isperfectpower`.
"""
function isprimepower(n::Integer)
    return isprime_power(n)
end

## Prime power function

#=
function __isprime(n)
	return isprime(fmpz(n))
end

function __factors(n)
	return n == 0 ? [] : factor(fmpz(n))
end
# function __factors(n)
# 	first(i) for i in Primes.factor(n)
# end

function primedivisors(n)
    n == 0 && return []
    __isprime(n) && return [fmpz(n)]
    f = __factors(n)
    return sort!([p for (p, e) ∈ f])
end

function ω(n)
	nprimedivisors = 0
	if n == 0
		nprimedivisors = 0
	elseif __isprime(n)
		nprimedivisors = 1
    else
        nprimedivisors = length(__factors(n))
	end

	return fmpz(nprimedivisors)

    # # original definition:
	# return fmpz(length(primedivisors(n)))
end

function Ω(n)
    n == fmpz(0) && return 0
    __isprime(n) && return fmpz(1)

    return sum(e for (__, e) ∈ __factors(n))
end

# this is a slow version compared to the Hecke one
function isperfectpower(n)
	return ω(n) == 1 && Ω(n) != 1
end

# this is slower than the Hecke verson
function isprimepower(n)
    return ω(n) == 1
end

=#

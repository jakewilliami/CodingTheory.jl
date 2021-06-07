#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

# TAKE ONE

using Nemo, Primes

function _isprime(n)
	return Nemo.isprime(Nemo.fmpz(n))
end

function _factors(n)
	return n == 0 ? [] : Nemo.factor(Nemo.fmpz(n))
end
# function _factors(n)
# 	first(i) for i in Primes.factor(n)
# end

function primedivisors(n)
    n == 0 && return []
    _isprime(n) && return [Nemo.fmpz(n)]
    f = _factors(n)
    return sort!([p for (p, e) ∈ f])
end

function ω(n)
	nprimedivisors = 0
	if n == 0
		nprimedivisors = 0
	end
	if _isprime(n)
		nprimedivisors = length(Nemo.fmpz(n))
	end

	return length(_factors(n))

	# return fmpz(length(primedivisors(n)))
end

function Ω(n)
    n == fmpz(0) && return 0
    _isprime(n) && return fmpz(1)
	
    return sum(e for (_, e) ∈ _factors(n))
end

function isperfectpower(n)
	return ω(n) == 1 && Ω(n) != 1
end

#---------------------------

#=
You find the smallest pair (a, n) (small with respect to n) such that a^n is your target number. Then you need to check if a is prime.
=#
function isprimepower(n::Integer)
	# return any(iszero(mod(n, p)) for p in primes(ceil(Int, sqrt(n))))
	
	n ≤ 1 && return false
	isprime(n) && return true
	ub = isqrt(n) + 1
		
	for a in primes(ub)
		for m in 2:ub # m > 1
			res = a^m
			res > n && break
			isequal(res, n) && return true
		end
	end
	
	return false
end

function primepower_brute_force(n::Integer)
	n ≤ 1 && return nothing
	isprime(n) && return n, 1
		
	for a in primes(ceil(Int, sqrt(n)))
		for m in 2:ceil(Int, sqrt(n)) # m > 1
			a^m > n && break
			isequal(a^m, n) && return a, m
		end
	end
	
	return nothing
end

function perfectpower_brute_force(n::Integer)
	n ≤ 1 && return nothing
	
	for a in 1:n
		for m in 2:n
			isequal(a^m, n) && return m, a
		end
	end
	
	return nothing
end

#========================================================================#

## The following function are obtained from DJB

struct DJBFloat
	a::Signed
	n::Unsigned
end

lg(x) = log2(x)

#=
From Knuths The Art of Computer Programming, semi-numerical algorithms, 2nd edition, 4.3.1, exercise 16.

Design an algorithm that divides a nonnegative n-place integer $(u_1u_2\ldots u_n)_b$ by $v$, where $v$ is a single-precision number (i.e., $0 < v < b$), producing the quotient $(w_1w_2\ldots 2_n)_b$ and remainder $r$.
=#

function Base.divrem(u::Integer, v::Real, b::Integer)
	r, j = 0, 1
	n = ndigits(u) #; base = b
	u̲ = round.(reverse(digits(u)), base = b)
	w̲ = Vector{Int}(undef, n)
	
	while j ≤ n
		w̲[j] = floor(Int, (r*b + u̲[j]) / v)
		r = mod(r*b + u̲[j], v)
		j += 1
	end
	
	# w = sum(w̲[i] * 10 ^ (length(w̲) - i) for i in eachindex(w̲))
	w = parse(Int, join(string.(w̲)))
	# w = foldr((i, j) -> 10 * j + i, reverse(w̲))
	#
	# n = 0
	# reduce(w̲) do i, j
	# 	10 * j + i
	# 	# global n = n + ndigits(w̲[i])
	# end
	# n = 0
	# w = sum(eachindex(w̲)) do i
	# 	w̲[i] * 10 ^ (length(w̲) - i + n)
	# 	global n = n + ndigits(w̲[i])
	# end
	# n = 0
	
	return w, r
end

Base.div(u::Integer, v::Real, b::Integer) = first(divrem(u, v, b))
Base.rem(u::Integer, v::Real, b::Integer) = last(divrem(u, v, b))

# function Base.div(a::Integer, n::Integer, k::Integer, b::Integer)
# 	f = nextpow(2, a)
# 	return a + f - ceil(Int, lg(k)) - b, floor(Int, n / (2^(f - ceil(Int, lg(k)) - b) * k))
# end

#=
In this section, I define truncation to $b$ bits, written $\text{trunc}_b$, and show that $\frac{r}{\text{trunc}_br}$ is between $1$ and $1+2^{1-b}$.  More generally, for any positive integer $k$, I define $\text{div}_b(r,k)$ as a floating point approximation to $\frac{r}{k}$, so that $\frac{r}{k\text{div}_b(r,k)}$ is between $1$ and $1+2^{1-b}$.\\ Fix $b\geq 1$.  Set $\text{div}_b(a,n,k)=(a+f-\left\lceil\text{lg} k\right\rceil - b, \left\lfloor\frac{n}{2^{f-\left\lceil\text{lg} k\right\rceil-b}k}\right\rfloor)$, where $2^{f-1}\leq n<w^f$.  (Note that $f-\left\lceil \text{lg} k\right\rceil - b$ may be negative.)  This map induces a map, also denoted $\text{div}_b$, upon positive floating-point numbers: \begin{equation}\text{div}_b(2^an,k)=2^{a+f-\left\lceil\text{lg} k\right\rceil - b}\left\lfloor\frac{n}{2^{f-\left\lceil\text{lg}k\right\rceil - b}k}\right\rfloor\quad\quad\quad\text{if }2^{f-1}\leq n<2^f.\end{equation}\\ To compute $\text{div}_b(r,k)$ I also use an algorithm designed for dividing by small integers; time spent computing $\text{div}_b$ is not $M$-time.
=#
function Base.trunc(n::Real, b::Integer)
    return n & (1 << b - 1)
	# return div(n, 1, b)
end


#=
Write a MIX program that multiplies $(u_1u_2\ldots u_n)_b$ by $v$, where $v$ is a single-precision number (i.e., $0 \leq v < b$), producing the answer $(w_0w_1\ldots w_n)_b$.  How much running time is required?

		ENTX	N		1	ADD		CARRY	N
		JOV		OFLO	1	JNOV	*+2		N
		ENTX	0		1	INCX	1		K
	2H	STX		CARRY	N	STA		W,1		N
		LDA		U,1		N	DEC1	1		N
		MUL		V		N	J1P		2B		N
		SLC		5		N	STX		W		1
		
The running time is 23N + K + 5 cycles, and K is rougly (1/2)N.
=#
function mul(u::Integer, v::Real, b::Integer)
	n = ndigits(u) #; base = b
	u̲ = round.(reverse(digits(u)), base = b)
	w̲ = Vector{Int}(undef, n)
end

#=
Given a positive floating-point number r and two positive integers b, k, to prime pow_b(r, k).  Let $\text{pow}_b(r,k)$ be the $b$-bit approximate $k$th power of $r$.
\begin{equation}\begin{aligned}
\text{pow}_b(r,1) &= \text{trunc}_br; \\
\text{pow}_b(r,2k) &= \text{trunc}_b(\text{pow}_b(r,k)^2); \\
\text{pow}_b(r,2k+1) &= \text{trunc}_b(\text{pow}_b(r,2k)\text{pow}_b(r,1))
\end{aligned}\end{equation}
=#
function __pow(r::Real, b::Integer, k::Integer)
	if r ≤ 0 || b ≤ 0 || k ≤ 0
		throw(error("r, b, and k must be positive numbers to compute pow_b(r, k)."))
	end
	println("pow: $r, $b, $k")
	isone(k) && return trunc(r, b)
	iszero(mod(k, 2)) && return trunc(__pow(r, b, k ÷ 2)^2, b) # use \div because k is even
	return trunc(__pow(r, b, k - 1) * trunc(r, b), b)
end

#=
I am trying to find a root $z$ of $z^ky-1$.
To compute nroot_b(y, k) for 1 \le b \le \ceil{lg8k}: In advance, find the exponent g satisfying 2^{g-  1} < y < 2^g, and set a = \floor{-g / k}, so that 2^a \le y^{-1/k} < 2^{a+1}.  Also set B = \ceil{lg(66(2k + 1))}
=#
function __nroot(y::Real, k::Real, b::Real)
	# find the exponent g satisfying 2^(g - 1) < y < 2^g
	g = nothing
	i = 0
	while true
		if 2^(float(i) - 1) < y && y ≤ float(2)^i
			g = i
			break
		end
		i += 1
	end
	
	a = floor(Int, -g / k)
	B = ceil(Int, lg(66 * (2*k + 1)))
	z = 2^float(a) + 2^(float(a) - 1)
	println("nroot: $z")
	j = 1
	
	while ! isequal(j, b)
		println("nroot: $y, $k, $j, $B")
		z = __nroot(y, k, j)
		r = trunc(__pow(z, k, B) * trunc(y, B), B)
		z = ifelse(r ≤ 993//1024, z + 2^(a - j + 1), z)
		z = ifelse(r > 1, z - 2^(a - j - 1), z)
		j = j + 1
	end
	
	return z
end

#=
To compute nroot_b(y, k) for b \ge 4 + \ceil{lg 8k} + 1: In advance set b' = \ceil{lg 2k} + \ceil{(b - \ceil{lg 2k}) / 2} and B = 2b' + 4 - \ceil{lg k}.  Note that b; < b.
=#
function __newton_nroot(y::Real, k::Real, b::Real)
	b′ = ceil(Int, (b - ceil(Int, lg(2*k))) / 2)
	B = 2*b′ + 4 - ceil(Int, lg(l))
	z = ifelse(b′ ≤ ceil(Int, lg(8*k)), __nroot(__nroot(y, k, b′)), __newton_nroot(y, k, b′)) # if else, then b′ ≥ ceil(Int, lg(8*k)) + 1
	r₂ = trunc(z, B) * (k + 1)
	r₃ = trunc(__pow(z, k + 1, B) * trunc(y, B), B)
	
	return r₄ = div(r₂ - r₃, k, B)
end

#=
Given positive integers n, x, k, to compute the sign of n - x^k: In advance, set f = \ceil{lg 2n}
=#
function __sign(n::Integer, x::Integer, k::Integer)
	if n ≤ 0 || x ≤ 0 || k ≤ 0
		throw(error("n, x, and k must be positive numbers to compute the sign of n - x^k using this algorithm."))
	end
	
	f = ceil(Int, lg(2n))
	b = 1
	while true
		r = __pow(x, k, b + ceil(Int, lg(8k)))
		
		n < r && return -1
		r*(1 + 2^(-b)) ≤ n && return 1
		b ≥ f && return 0
		
		b = min(2*b, f)
	end
end

#=
Given integers n \ge 2 and k \ge 2, and a positive floating-point number y, to see if n is a kth power: in advance, set f = \floor{lg 2n} and b = 3 + \ceil{f / k}
=#
function __iskthpower(n::Integer, k::Integer, y::Real)
	if n < 2 || k < 2 || y ≤ 0
		thow(error("n, k, and y must be positive (and n and k must be more than one) to see if n is a kth power using this algorithm."))
	end
	
	f = floor(Int, lg(2*n))
	b = 3 + ceil(Int, f / k)
	r = __nroot(y, k, b)
	
	x = nothing
	for i in 1:r
		if abs(r - x) ≤ 5//8
			x = i
		end
	end
	
	(iszero(x) || abs(r - x) ≥ 1//4) && return 0
	
	sign = __sign(n, x, k) # why do we do this step?!
	isequal(n, x^k) && return x
	
	return 0
end

#=
Given an integer n \ge 2, to decompose n as a perfect power if possible: in advance, set f = \floor{lg 2n}
=#
function __find_prime_power(n::Integer)
	if n < 2
		throw(error("This algorithm to compute the prime power if possible only works for integers more than one."))
	end
	
	# isequal(n, 2) && return n, 1
	
	f = floor(Int, lg(2*n))
	y = __nroot(n, 1, 3 + ceil(Int, f / 2))
	
	for p in primes(f)
		x = __iskthpower(n, p, y)
		x > 0 && return x, p
	end
	
	return x, 1
end

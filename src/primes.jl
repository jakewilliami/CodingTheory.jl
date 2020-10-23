#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
using Primes: primes, isprime

#=
You find the smallest pair (a,n) (small with respect to n) such that a^n is your target number. Then you need to check if a is prime.
=#
function isprimepower(n::Integer)::Bool
	isprime(n) && return true
	
	for a in primes(Int(floor(sqrt(n))))
		for m in 1:Int(floor(sqrt(n)))
			isequal(a^m, n) && return true
		end
	end
	
	return false
end

function isprimepower_log(n::Integer)::Bool
	for p in prime(Int(floor(sqrt(n))))
		if isinteger(log(n, p))
			return true
		end
	end
	
	return false
end

## The following function are obtained from https://cr.yp.to/papers/powers.pdf

#=
Log base 2: alternative notation
=#
lg(x::Number) = log2(x)

#=
In this section, I define truncation to $b$ bits, written $\text{trunc}_b$, and show that $\frac{r}{\text{trunc}_br}$ is between $1$ and $1+2^{1-b}$.  More generally, for any positive integer $k$, I define $\text{div}_b(r,k)$ as a floating point approximation to $\frac{r}{k}$, so that $\frac{r}{k\text{div}_b(r,k)}$ is between $1$ and $1+2^{1-b}$.\\ Fix $b\geq 1$.  Set $\text{div}_b(a,n,k)=(a+f-\left\lceil\text{lg} k\right\rceil - b, \left\lfloor\frac{n}{2^{f-\left\lceil\text{lg} k\right\rceil-b}k}\right\rfloor)$, where $2^{f-1}\leq n<w^f$.  (Note that $f-\left\lceil \ltext{lg} k\right\rceil - b$ may be negative.)  This map induces a map, also denoted $\text{div}_b$, upon positive floating-point numbers: \begin{equation}\text{div}_b(2^an,k)=2^{a+f-\left\lceil\text{lg} k\right\rceil - b}\left\lfloor\frac{n}{2^{f-\left\lceil\text{lg}k\right\rceil - b}k}\right\rfloot&\text{if }2^{f-1}\leq n<2^f.\end{equation}\\ To compute $\text{div}_b(r,k)$ I also use an algorithm designed for dividing by small integers; time spent computing $\text{div}_b$ is not $M$-time.
=#
function Base.div(a::Integer, n::Integer, k::Integer, b::Integer)
	a + f - ceil(lg(k)) - b, n / (2^(f - ceil(lg(k)) - b) * k)
end


#=
Given a positive floating-point number r and two positive integers b, k, to prime pow_b(r, k).
=#
function __pow(r::Real, b::Integer, k::Integer)
	if r ≤ 0 || b ≤ 0 || k ≤ 0
		throw(error("r, b, and k must be positive numbers to compute pow_b(r, k)"))
	end
	
	isone(k) && return trunc(r, base = b)
	
	if iszero(mod(k, 2))
		__pow(r, b, k / 2)
		return trunc(__pow(r, b, k / 2)^2, base = b)
	end
	
	__pow(r, b, k - 1)
	return trunc(__pow(r, b, k - 1) * trunc(r, base = b), base = b)
end

#=
To compute nroot_b(y, k) for 1 \le b \le \ceil{lg8k}: In advance, find the exponent g satisfying 2^{g-  1} < y < 2^g, and set a = \floor{-g / k}, so that 2^a \le y^{-1/k} < 2^{a+1}.  Also set B = \ceil{lg(66(2k + 1))}
=#
function __nroot(y::Real, k::Real, b::Real)
	# find the exponent g satisfying 2^(g - 1) < y < 2^g
	g = nothing
	for i in 1:y
		if 2^(i - 1) < y && y < 2^i
			g = i
			break
		end
	end
	
	a = floor(-g / k)
	B = ceil(lg(66 * (2*k + 1)))
	z = 2^a + 2^(a - 1)
	j = 1
	
	nroot(y, k, j) = z
	while ! isequal(j, b)
		r = trunc(__pow(z, k, B) * trunc(y, base = B), base = B)
		
		if r ≤ 993//1024	z = z + 2^(a - j + 1)	end
		if r > 1	z = z - 2^(a - j - 1)	end
		
		j = j + 1
	end
end

#=
To compute nroot_b(y, k) for b \ge 4 + \ceil{lg 8k} + 1: In advance set b' = \ceil{lg 2k} + \ceil{(b - \ceil{lg 2k}) / 2} and B = 2b' + 4 - \ceil{lg k}.  Note that b; < b.
=#
function __newton_nroot(y::Real, k::Real, b::Real)
	b′ = ceil((b - ceil(lg(2*k))) / 2)
	B = 2*b′ + 4 - ceil(lg(l))
	z = nothing
	
	if b′ ≤ ceil(lg(8*k))
		z = __nroot(__nroot(y, k, b′))
	elseif b′ ≥ ceil(lg(8*k)) + 1 # else
		z = __newton_nroot(y, k, b′)
	end
	
	r2 = trunc(z, base = B) * (k + 1)
	r3 = trunc(__pow(z, k + 1, B) * trunc(y, base = B), base = B)
	
	return r4 = div(r2 - r3, k, base = B)
end

#=
Given positive integers n, x, k, to compute the sign of n - x^k: In advance, set f = \ceil{lg 2n}
=#
function __sign(n::Integer, x::Integer, k::Integer)
	if n ≤ 0 || x ≤ 0 || k ≤ 0
		throw(error("n, x, and k must be positive numbers to compute the sign of n - x^k using this algorithm."))
	end
	
	f = ceil(lg(2n))
	b = 1
	while true
		r = __pow(x, k, b + ceil(lg(8k)))
		
		n < r && return -1
		r*(1 + 2^(-b)) ≤ n && return 1
		b ≥ f && return 0
		
		b = min(2*b, f)
	end
end

#=
Given integers n \ge 2 and k \ge 2, and a positive floating-point number y, to see if n is a kth power: in advance, set f = \floor{lg 2n} and b = 3 + \ceil{f / k}
=#
function iskthpower(n::Integer, k::Integer, y::Real)
	if n < 2 || k < 2 || y ≤ 0
		thow(error("n, k, and y must be positive (and n and k must be more than one) to see if n is a kth power using this algorithm."))
	end
	
	f = floor(lg(2*n))
	b = 3 + ceil(f / k)
	r = __nroot(y, k, b)
	
	x = nothing
	for i in 1:r
		if abs(r - x) ≤ 5//8
			x = i
		end
	end
	
	if	iszero(x) || abs(r - x) ≥ 1//4	return 0	end
	
	sign = __sign(n, x, k)
	isequal(n, x^k) && return x
	
	return 0
end

#=
Given an integer n \ge 2, to decompose n as a perfect power if possible: in advance, set f = \floor{lg 2n}
=#
function __find_prime_power(n::Integer)
	if n ≤ 1
		throw(error("This algorithm to compute the prime power if possible only works for integers more than one."))
	end
	
	f = floor(lg(2*n))
	y = __nroot(n, 1, 3 + ceil(f / 2))
	
	for p in primes(f)
		x = iskthpower(n, p, y)
		x > 0 && return x, p
	end
	
	return (x, 1)
end

### TAKE TWO

using Primes: primes, factor, isprime

### Helpers

nroot(x::Number, n::Number) = x^(1/n)

function factors(n::Integer)
    f = [one(n)]

    for (p,e) in factor(n)
        f = reduce(vcat, [f*p^j for j in 1:e], init=f)
    end

    return ifelse(isone(length(f)), [one(n), n], sort!(f))
end
 
"""
    ord(q::Integer, n::Integer) -> Integer

Computes the multaplicative order, ord_q n, of a number n mod q.
"""
function ord(q::Integer, n::Integer)
    isone(gcd(n, q)) || throw(error("$n and $q are not coprime"))
    res = one(q)

    for (p, e) in factor(q)
        m = p^e
        t = div(q, p) * (p-1)
        for f in factors(t)
            if isone(powermod(n, f, q))
                res = lcm(res, f)
                break
            end
        end
    end

    return res
end



### Main code
#=
\begin{enumerate} \item For each prime power $q$ such that $2^q\leq n$, write down a positive integer $r_q$ such that if $n$ is a $q$\th power, then $n=r^q_q$.\item Find a finite coprime set $P$ of integers larger than 1 such that each of $n,r_2,r_3,r_4, r_5,r_7,\ldots$ is a product of powers of elements of $P$.  (In this paper, ``coprime'' means ``pairwise coprime.'')\item Factor $n$ as $\Pi_{p\in P}p^{n_p}$, and compute $k=\text{gcd}n_p:p\in P$. \end{enumerate}=#

#=
From \texttt{Powers2.pdf}:\begin{quotation}
	Define $r_q$ as an integer 2-adically close to $n^{1/q}$, as explained in \texttt{Powers.pdf}, section 21.\end{quotation}

From \texttt{Powers.pdf}, section 21:\begin{quotation}
	In this section, I describe a \emph{2-adic variant} of Algorithm X.  With this variant, I can work with integers rather than floating point numbers; I no longer need to guard bits; I can jump directly into Newton's method without a preliminary binary search; and a proper error analysis takes a few lines rather than several pages.\par
	It will be convenient to restrict attention to odd $n$.  See section 22 to handle even $n$.\par
	\textbf{Motivation.}  To check if $n$ is a $k$\th power, I compute a tentative $k$\th root of $n$\emdash an integer x such that no integer other than $x$ can possibly be the $k$\th root of $n$.  Then I test whether $x^k=n$.\par
	To find $x$ in Part II and Part III, I constructed a number that was \emph{close} to a $k$\th root of $n$ in the usual metric.  I used the same metric again to check whether $x^k=n$: I computed $x^k$ in low precision to see whether it was close to $n$.\par
	Nothing in the original problem suggests this metric.  The 2-adic variant uses a different metric, where $i$ and $j$ are close if $i-j$ is divisible by a high power of 2.\par
	\textbf{Notation.}  This section deviates from the notation of parts I, II, and III: $r$, $y$, and $z$ are odd integers rather than positive floating-point numbers.\par
	\textbf{Lemma 21.1.}  \emph{If $2i\equiv 2j \mod 2^{b+1}$, and $b\geq 1$, then $i^2\equiv j^2\mod w^{b+1}.$}\par
	\textbf{2-adoc approximate powers.}  Fix positive integers $k$ and $b$.  For any integer $m$ define $\text{pow}_{2,b}(m,k)=m^k\mod 2^b$.  See section 6 for methods of computing $\text{pow}_{2, b}$ without many multiplications.  As I compute $\text{pow}_{2,b}(m, k)$, I keep track of an ``overflow bit'' to figure out whether $m^k\mod 2^b=m^k$.
\end{quotation}
	
From \texttt{Powers.pdf}, section 6:\begin{quotation}
	Let $r$ be a positive floating-point number, and let $k$ be positive integers.  Then $\text{pow}_{b}(r, k)$, the \textbf{$\mathbf{b}$-bit approximate $\mathbf{k}th$ power of $\mathbf{r}$} is a floating-point approximation to $r^k$.  In this section, I show how to compute $\text{pow}_b(r,k)$ in $M$-time at most $P(k)M(b), where $P(k)\leq 2\lfloor\text{lg} k\rfloor$.\par
	Define $P(k)$ for $k\geq 1$ as follows: $P(1)=0;P(2k)=P(k)+1;P(2k+1)=P(2k)+1$.\par
	\textbf{Lemma 6.1.}  $P(k)\leq 2\lfloor\text{lg} k\rfloor$.\par
	\emph{Proof.} $\ldots$\par
	\textbf{Lemma 6.2} $\text{pow}_b(r,k)\leq r^k<\text{pow}_b(r,k)(1+2^{1-b})^{2k-1}.$\par
	\emph{Proof.} $\ldots$\par
	\emph{Algorithm P.} $\ldots$\par
	\textbf{Lemma 6.3.}  Algorithm P computes $\text{pow}_b(r, k)$ in $M$-time at most $P(k)M(b)$.\par
	\emph{Proof.} $\ldots$\par
\end{quotation}

From \texttt{Powers.pdf}, section 22:\begin{quotation}
	
\end{quotation}
=#


lg(x) = log2(x)

#=
From Knuths The Art of Computer Programming, semi-numerical algorithms, 2nd edition, 4.3.1, exercise 16.

Design an algorithm that divides a nonnegative n-place integer $(u_1u_2\ldots u_n)_b$ by $v$, where $v$ is a single-precision number (i.e., $0 < v < b$), producing the quotient $(w_1w_2\ldots 2_n)_b$ and remainder $r$.
=#

"""
	divrem(b::Integer, u::Integer, v::Real)

This algorithm divides a non-negative n-place integer (u1 u2 ... un) by v, where v is a single-precision number (i.e., 0 < v < b), producing the quotient (w1 w2 ... wn) and remainder r.

This algorithm was written by Donald Knuth, in The Art of Computer Programming, semi-numerical algorithms, 2nd edition, exercise 4.3.1, 16.
"""
function Base.divrem(b::Integer, u::Integer, v::Real)
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

"""
	div(b::Integer, u::Integer, v::Real)

This algorithm divides a non-negative n-place integer (u1 u2 ... un) by v, where v is a single-precision number (i.e., 0 < v < b), producing the quotient (w1 w2 ... wn) and ommiting the remained.  To return the quotient and the remainder, see `divrem`.
"""
Base.div(b::Integer, u::Integer, v::Real) = first(divrem(b, u, v))

"""
	rem(b::Integer, u::Integer, v::Real)

This algorithm divides a non-negative n-place integer (u1 u2 ... un) by v, where v is a single-precision number (i.e., 0 < v < b), producing a remainder r.  To return the quotient and the remainder, see `divrem`.
"""
Base.rem(b::Integer, u::Integer, v::Real) = last(divrem(b, u, v))

#=
In this section, I define truncation to $b$ bits, written $\text{trunc}_b$, and show that $\frac{r}{\text{trunc}_br}$ is between $1$ and $1+2^{1-b}$.  More generally, for any positive integer $k$, I define $\text{div}_b(r,k)$ as a floating point approximation to $\frac{r}{k}$, so that $\frac{r}{k\text{div}_b(r,k)}$ is between $1$ and $1+2^{1-b}$.\\ Fix $b\geq 1$.  Set $\text{div}_b(a,n,k)=(a+f-\left\lceil\text{lg} k\right\rceil - b, \left\lfloor\frac{n}{2^{f-\left\lceil\text{lg} k\right\rceil-b}k}\right\rfloor)$, where $2^{f-1}\leq n<w^f$.  (Note that $f-\left\lceil \text{lg} k\right\rceil - b$ may be negative.)  This map induces a map, also denoted $\text{div}_b$, upon positive floating-point numbers: \begin{equation}\text{div}_b(2^an,k)=2^{a+f-\left\lceil\text{lg} k\right\rceil - b}\left\lfloor\frac{n}{2^{f-\left\lceil\text{lg}k\right\rceil - b}k}\right\rfloor\quad\quad\quad\text{if }2^{f-1}\leq n<2^f.\end{equation}\\ To compute $\text{div}_b(r,k)$ I also use an algorithm designed for dividing by small integers; time spent computing $\text{div}_b$ is not $M$-time.
=#

"""
	trunc(b::Integer, n::Real)
	
Truncates a number `n` to `b` bits.
"""
function Base.trunc(b::Integer, n::Real)
	return n & (1 << b - 1)
	# return div(b, n, 1)
end

#=
Given a positive floating-point number r and two positive integers b, k, to prime pow_b(r, k).  Let $\text{pow}_b(r,k)$ be the $b$-bit approximate $k$th power of $r$.
\begin{equation}\begin{aligned}
\text{pow}_b(r,1) &= \text{trunc}_br; \\
\text{pow}_b(r,2k) &= \text{trunc}_b(\text{pow}_b(r,k)^2); \\
\text{pow}_b(r,2k+1) &= \text{trunc}_b(\text{pow}_b(r,2k)\text{pow}_b(r,1))
\end{aligned}\end{equation}
=#

"""
	pow(b::Integer, r::Real, k::Integer)
	
Given a positive floating-point number r and two positive integers b and k, pow(b, r, k) computes the b-bit approximate kth power of r.  pow_b(r, k) ≤ r^k < pow_b(r, k)(1 + 2^(1 - b))^(2k - 1).
"""
function pow(b::Integer, r::Real, k::Integer)
	if r ≤ 0 || b ≤ 0 || k ≤ 0
		throw(error("r, b, and k must be positive numbers to compute pow_b(r, k)."))
	end
	isone(k) && return trunc(b, r)
	iszero(mod(k, 2)) && return trunc(b, pow(b, r, k ÷ 2)^2) # use \div because k is even
	return trunc(b, pow(b, r, k - 1) * trunc(b, r))
end

"""
	pow2(b::Integer, r::Real, k::Integer)
	
Given a positive floating-point number r and two positive integers b and k, pow2(b, r, k) computes the 2-adic b-bit approximate kth power of r.
"""
function pow2(b::Integer, m::Integer, k::Integer)
	if k ≤ 0 || b ≤ 0
		throw(error("k, and b must be positive integers to find a 2-adic approximate power."))
	end
		
	return mod(m^k, 2^b)
end

"""
	tentative_check_kth_roots(n::Integer, x::Integer, k::Integer)

A straight forward algorithm for checking whether x^k = n.
"""
function tentative_check_kth_roots(n::Integer, x::Integer, k::Integer)
	if n ≤ 0 || x ≤ 0 || k ≤ 0
		throw(error("n, x, and k must be positive integers to compute check whether n=x^k."))
	end
	
	f = floor(Integer, lg(2*n))
	
	isone(x) && return ifelse(isone(n), true, false) # return ifelse(isone(n), 0, 2)
	
	b = 1
	while true
		# r = pow(x, b, k) # 2, b
		r = pow2(b, x, k)
		isequal(mod(n, 2^b), r) || return false # 2
		b ≥ f && return ifelse(isequal(r, x^k), true, false) # ifelse(isequal(r, x^k), 0, 2)
		
		b = min(2*b, f)
	end
end

#=
\textbf{Lemma 21.2}.  Algorithm C2 prints 0 if and only if n=x^k.\par
\emph{Proof.} $\ldots$\par
\textbf{2-adic approximate multiplication and division}.  Fix $b\geq 1$.  For $m$ an integer and $k$ a positive integer write $\text{mul}_{2, b}(m, k) = km \mod 2^b.$\par
If $k$ is odd, write $\text{div}_{2, b}(m, k)$ for the unique integer between 0 inclusive and $2^b$ exclusive such that $m\equiv k \text{div}_{2, b}(m, k) \mod 2^b$.
=#

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

"""
	mul2(b::Integer, m::Integer, k::Integer)

The notation mul(r, k) is used in the Bernstein paper to denote k * r.  The notation mul_{2, b}(m, k) represents 2-adic multiplication.
"""
function mul2(b::Integer, m::Integer, k::Integer)
	if k ≤ 0
		throw(error("Cannot compute mul(m, k) when k is negative using the 2-adic method."))
	end
	
	return mod(k*m, 2^b)
end

"""
	div2(b::Integer, m::Integer, k::Integer)
div_b is defined as a floating-point approximation to r/k, so that r/k div_b(r, k) is between 1 and 1 + 2^(1 - b).  div2 is the 2-adic approximation of this.
"""
function div2(b::Integer, m::Integer, k::Integer)
	if k ≤ 0
		throw(error("Cannot compute mul(m, k) when k is negative using the 2-adic method."))
	end
	
	for m in 0:2^b
		if isequal(m, mod(div(b, m, k), 2^b))
			return div(b, m, k)
		end
	end
end

"""
	odd_nroot2(b::Integer, y::Integer, k::Integer)

Fix an odd integer y and a positive odd integer k.  odd_nroot2 finds an approximate negative kth root of y by Newton's method.
"""
function odd_nroot2(b::Integer, y::Integer, k::Integer)
	if k ≤ 0 || b ≤ 0
		throw(error("Cannot use this method for negative number k and b."))
	end
	
	if iseven(y) || iseven(k)
		throw(error("Cannot use this method for even y and k."))
	end
	
	b′ = ceil(Integer, b / 2)
	
	isone(b) && return 1
	z = odd_nroot2(b′, y, k)
	r₂ = mul2(b, z, k + 1)
	r₃ = mod(y * pow2(b, z, k + 1), 2^b)
	return r₄ = div2(b, r₂ - r₃, k)
end

#=
\textbf{Lemma 21.3.}  If $k$ is odd and $r=\text{nroot}_{2,b}(y, k)$, then $r^k y\mod 2^b = 1$.\par
\emph{Proof.} $\ldots$\par
=#

"""
	sqrt_nroot2(b::Integer, y::Integer)

Fix odd integer y.  For each b ≥ 1, define and construct nroot_(2, b)(y, 2).
"""
function sqrt_nroot2(b::Integer, y::Integer)
	# if b ≤ 0
	# 	throw(error("Cannot use this method for negative number b."))
	# end
	#
	# if iseven(y)
	# 	throw(error("Cannot use this method for even y."))
	# end
	
	b′ = ceil(Integer, (b + 1) / 2)
	
	isone(b) && return ifelse(isone(mod(y, 4)), 1, 0)
	isequal(b, 2) && return ifelse(isone(mod(y, 8)), 1, 0)
	iszero(z) && return 0
	
	r₂ = mul2(b + 1, z, 3)
	r₃ = mod(y * pow2(b + 1, z, 3), 2^(b + 1))
	return r₄ = mod(r₂ - r₃, 2^b)
end

#=
\textbf{Lemma 21.4.}  Set $r = \text{nroot}_{2, b}(y, 2)$.  If $i^2y\mod 2^{b + 1}=1$ for some odd integer $i$, then $r\neq 0$.  If $r\neq 0$ then $r^2y\mod 2^{b+1} = 1$.\par
\emph{Proof.} $\ldots$\par
=#

"""
	perfect_power_decomp_odd(n::Integer, k::Integer, y::Integer)

Given a positive odd integer n, an integer k ≥ 2 such that either k = 2 or k is odd, and an odd integer y, `perfect_power_decomp_odd` checks if n is a kth power.
"""
function perfect_power_decomp_odd(n::Integer, k::Integer, y::Integer)
	if iseven(n) && !isequal(k, 2)
		throw(error("Cannot use this algorithm on even n."))
	end
	
	if k < 2
		throw(error("k must be ≥ 2."))
	else
		if iseven(k) && !isequal(k, 2)
			throw(error("k must be 2 or odd."))
		end
	end
	
	f = floor(Integer, lg(2*n))
	b = ceil(Integer, f / k)
	
	r = odd_nroot2(b, y, k)
	(isequal(k, 2) && iszero(r)) && return 0
	tentative_check_kth_roots(n, r, k) && return r
	(isequal(k, 2) && tentative_check_kth_roots(n, 2^b - r, k)) && return 2^b - r
	return 0
end

#=
\textbf{Lemma 21.5.}  Set $f=\lfloor\text{lg}2n\rfloor$ and $b=\lceilf/k\rceil$.  Assume that $yn\mod 2^{b+1}=1$.  If $n$ is a $k$\th power, Algorithm K2 (above) prints $n^{1/k}$.  If $n$ is not a $k$\th power, Algorithm K2 prints 0.\par
\emph{Proof.} $\ldots$\par
=#

"""
	odd_decomp_perf_pow(n::Integer)

Given an odd integer n ≥ 2, `odd_decomp_perf_pow` decomposes n as a perfect power if possible.
"""
function perfectpower(n::Integer)
    # isprime(n) && return n, 1

	# if iseven(n) && !isequal(n, 2)
	# 	throw(error("This method must have odd n."))
	# end
	
	if n < 2
		throw(error("n must be ≥ 2."))
	end
    
    isprime(n) && return n, 1

    # sqrt_val = sqrt(n)
    # isinteger(sqrt_val) && return round(Integer, sqrt_val), round(Integer, sqrt_val)

    if iseven(n)
        if ispow2(n)
            return 2, round(Integer, log2(n))
        else
            return n, 1
        end
    end

	f = floor(Integer, lg(2*n))
	
	y = odd_nroot2(ceil(Integer, f / 2) + 1, n, 1)

	for p in primes(f - 1)
		x = perfect_power_decomp_odd(n, p, y)
		x > 0 && return x, p
	end
	
	return n, 1
end

#=
\textbf{Lemma 21.6.}  If $n$ is a perfect power, the above algorithm (Algorithm X2) prints a prime number $p$ and a positive integer $x$ such that $x^p=n$.  If $n$ is not a perfect power, Algorithm X2 print (n, 1).
\emph{Proof.} $\ldots$\par
=#

#=
As usual, fix $n\geq 2$.  In this section I discuss several tricks based on computing $n$ mod $q$ for on or more primes $q$.\par
\textbf{If $n$ has no small prime divisors, lower the exponent bound.}  If $n$ is odd and $n=x^k$, then $x$ is also odd.  So $x\geq 3$ and $k\leq \log_3 n$.  More generally, one may compute $n\mod q$ for all primes $p < T$, for some bound $T$.  If $n\mod q$ is always nonzero, one need not check exponents past $\log_t n$.\par
\textbf{If n has a prime divisor, find its order.}  What if $n\mmod 1 = 0$?  First compute the number $\text{ord}_q n$ of factors $q$ in $n$, together with $n/q^{\text{ord}_qn}$.  Then check, for each $p$ dividing $\text{ord}_qn$ whether $n/q^{\text{ord}_qn}$ is a $p$\th power.  Otherwise, $n$ cannot be a perfect power.  (Note that $n/q^{\text{ord}_qn}$ may be 1, in which case no trsting is necessary.)\par Recall that the 2-adic method in section 21 requires $n$ to be odd.  This is not a serious restriction.  If $n$ is even, the method here ends up checking whether $n/2^{\text{ord}_2n}$ is a $p$th power, for various primes $p$; and $n/2^{\text{ord}_2n}$ is odd.\par There are several plausible ways to compute the number of factors $q$ in $n$.\par If $q=2$ then $\text{ord}_qn$ is the number of 0 buts at the bottom of $n$'s binary expansion.  \par If $q>2$, I do a binary search upon the number of factors.  The idea is to compute $n\mod q^c$ and $\lfloor n/q^c\rfloor$ for some integers $c\approx (\log_q n)/2$.  If $n\mod q^c\neq 0$ then $\text{ord}_qn=\text{ord}_q(n\mod q^c)$; if $n\mod q^c=0$ then $\text{ord}_qn=c+\text{ord}_q(n/q^c)$.  Chop $c$ in half and repeat.  This method takes essentially linear time given fast multiplication and the algorithms from section 8.\par I could instead do a linear search; this amounts to always taking $c=1$ in the above description.  This will be faster than a binary search on average.  I could compromise with a sequence of $c$ that is at first optimistic but backs off quickly if necessary.  For example: begin with $c=1$, double $c$ if $n\mod q^c = 0$, and chop $c$ in half if $n\mod q^c\neq 0$.\par
\textbf{Check the character of residues of n}.  If $n$ is a $k$\th power, and $q$ is a prime with $q\mod k=1$, then $n^{(q - 1) / k}\mod q$ is either 0 or 1.  A non-$k$\th power has roughly a $1/k$ chance of passing this test.\par
\textbf{Check the residues of tentative roots.}  If $n=x^k$ then $n\mod q=x^k\mod q$.  So, given $n\mod q$, I can try to weed out a tentative root $x$ by calculating the $k\th$ power modulo $q$ of $x\mod q$.  In practice this test is quite powerful: fi $n\neq x^k$ then very few primes $q$ divide $n-x^k$.\par In this test $q$ need not be prime.  It might be convenient to check whether $n$ agrees with $x^k$ modulo $2^{32}$, for example, though this is redundant if $x$ was constructed by 2-adic methods.\par One could develop a fast randomized power-testing algorithm along these lines.  Start from a tentative root $x$.  First check if $x^k\leq 2n$.  Then check if $n\mod q$ equals $x^k\mod q$ for a set of ``random'' primes $q$ with product larger than $n$.  This test will succeed if and only if $n=x^k$.  If $n\neq x^k$ then one will, on average, test very few $q$'s.\par
\textbf{Check for small divisors of tentative roots}.  If $n$ is not divisible by any primes $q<T$, and $n=x^k$, then $x$ is not divisible by any primes $q<T$.  So I can throw away any tentative root $x$ that has prime factors smaller than $T$.  This is much weaker than testing whether $n\mod d=x^k\mod q$ for each $q<T$, but it is also much faster.
=#

function lower_exp_bound(n::Integer, x::Integer, k::Integer)
    # if isodd(n)
        # continue
    # end
end

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

function isperfectpower(n)
	return ω(n) == 1 && Ω(n) != 1
end

function isprimepower(n)
    return ω(n) == 1
end

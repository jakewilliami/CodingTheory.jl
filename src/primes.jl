#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#
    
using Primes

function isprimepower(n::Integer)::Bool
    if isone(n)
        return false
    end
    
    if isequal(n, 2) || isequal(n, 3) || isequal(n, 4)
        return true
    end
    
    if Primes.isprime(n)
        return true
    end
    
    for i in Primes.primes(2, Int(ceil(sqrt(n))))
        m = n / i
        if isinteger(m)
            n = m
        end
        
        if isone(n)
            return true
        end
    end
    
    return false
end

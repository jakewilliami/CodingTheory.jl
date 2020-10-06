#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

include(joinpath(dirname(dirname(@__FILE__)), "src", "CodingTheory.jl"))

using .CodingTheory
using Test
using Polynomials

@testset "CodingTheory.jl" begin
	p = Polynomial([1, 1, 2, 0, 1, 2, 1])
	q = Polynomial([2, 1, 1])
	a = Polynomial([0, 1, 0, 1, 0, 0, 1])
	b = Polynomial([1, 0, 1])

	@test mod(rem(p, q), 3) == Polynomial([1, 1])
	@test mod(rem(a, b), 2) == Polynomial([1])
	
	@test multiplication_table(2, 3) == Polynomial[Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]); Polynomial([0]) Polynomial([1]) Polynomial([2]) Polynomial([0, 1]) Polynomial([1, 1]) Polynomial([2, 1]) Polynomial([0, 2]) Polynomial([1, 2]) Polynomial([2, 2]); Polynomial([0]) Polynomial([2]) Polynomial([1]) Polynomial([0, 2]) Polynomial([2, 2]) Polynomial([1, 2]) Polynomial([0, 1]) Polynomial([2, 1]) Polynomial([1, 1]); Polynomial([0]) Polynomial([0, 1]) Polynomial([0, 2]) Polynomial([0, 0, 1]) Polynomial([0, 1, 1]) Polynomial([0, 2, 1]) Polynomial([0, 0, 2]) Polynomial([0, 1, 2]) Polynomial([0, 2, 2]); Polynomial([0]) Polynomial([1, 1]) Polynomial([2, 2]) Polynomial([0, 1, 1]) Polynomial([1, 2, 1]) Polynomial([2, 0, 1]) Polynomial([0, 2, 2]) Polynomial([1, 0, 2]) Polynomial([2, 1, 2]); Polynomial([0]) Polynomial([2, 1]) Polynomial([1, 2]) Polynomial([0, 2, 1]) Polynomial([2, 0, 1]) Polynomial([1, 1, 1]) Polynomial([0, 1, 2]) Polynomial([2, 2, 2]) Polynomial([1, 0, 2]); Polynomial([0]) Polynomial([0, 2]) Polynomial([0, 1]) Polynomial([0, 0, 2]) Polynomial([0, 2, 2]) Polynomial([0, 1, 2]) Polynomial([0, 0, 1]) Polynomial([0, 2, 1]) Polynomial([0, 1, 1]); Polynomial([0]) Polynomial([1, 2]) Polynomial([2, 1]) Polynomial([0, 1, 2]) Polynomial([1, 0, 2]) Polynomial([2, 2, 2]) Polynomial([0, 2, 1]) Polynomial([1, 1, 1]) Polynomial([2, 0, 1]); Polynomial([0]) Polynomial([2, 2]) Polynomial([1, 1]) Polynomial([0, 2, 2]) Polynomial([2, 1, 2]) Polynomial([1, 0, 2]) Polynomial([0, 1, 1]) Polynomial([2, 0, 1]) Polynomial([1, 2, 1])]
	
	@test list_span([2, 1, 1], [1, 1, 1], 3) == [[0, 0, 0], [1, 1, 1], [2, 2, 2], [2, 1, 1], [0, 2, 2], [1, 0, 0], [1, 2, 2], [2, 0, 0], [0, 1, 1]]
	
	@test islinear([[0,0,0],[1,1,1],[1,0,1],[1,1,0]], 2) == false
	@test islinear([[0,0,0],[1,1,1],[1,0,1],[0,1,0]], 2) == true
	
	@test code_distance([[0,0,0,0,0],[1,0,1,0,1],[0,1,0,1,0],[1,1,1,1,1]]) == 2
	@test code_distance([[0,0,0,0,0],[1,1,1,0,0],[0,0,0,1,1],[1,1,1,1,1],[1,0,0,1,1],[0,1,1,0,0]]) == 1
	
	@test Alphabet("123").Σ == [1, 2, 3]
	@test Alphabet([1, 2, 3]).Σ == [1, 2, 3]
	@test Alphabet(["1", "2", "3"]).Σ == [1, 2, 3]
	
end # end runtests

# @btime test()

# C1 = [[0,0,0,0,0],[1,0,1,0,1],[1,0,0,0,1],[0,1,0,1,0]]
# C2 = [[0,0,0,0,0],[1,1,1,0,0],[0,0,0,1,1],[1,1,1,1,1],[1,0,0,1,1],[0,1,1,0,0]]
# C3 = [[0,0,0,0,0],[1,0,1,0,1],[0,1,0,1,0],[1,1,1,1,1]]

# println(islinear(C1, 2))#addition
# println(islinear(C2, 2))#addition
# println(islinear(C3, 2))

# println(code_distance(C1))
# println(code_distance(C2))
# println(code_distance(C3))

# @show list_span([1, 0, 1, 0, 1, 0], [0, 1, 0, 0, 1, 0], [1, 1, 1, 1, 1, 1], 2)

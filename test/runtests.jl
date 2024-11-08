using CodingTheory
import CodingTheory.deepsym, CodingTheory.displaymatrix

using Test
using Polynomials

function aresamearrays(A::AbstractVector{T}, B::AbstractVector{R}) where {T, R}
    A === B && return true
    length(A) ≠ length(B) && return false

    for (a, b) in zip(sort(A), sort(B))
        if a ≠ b
            return false
        end
    end

    return true
end

"""
Given two arrays A and B, ensures they have the same elements
"""
function aresamearrays(A::AbstractArray{T, N}, B::AbstractArray{R, M}) where {T, R, N, M}
    A === B && return true
    length(A) ≠ length(B) && return false

    for (a, b) in zip(sort(reshape(A, :)), sort(reshape(B, :)))
        # for (a, b) in zip(sort(A, dims = 1), sort(B, dims = 1))
        if a ≠ b
            return false
        end
    end

    return true
end

@time @testset "CodingTheory.jl" begin
    @test has_identity([1 0 0 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0], 1) == true
    @test has_identity([1 0 0 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]) == true
    @test has_identity([1 0 1 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]) == true
    @test has_identity([1 0 1 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0], 3) == false
    @test has_identity([-96 -66 20 1 0 0; -65 59 -82 0 1 0; -16 87 -113 0 0 1]) == true
    @test has_identity([
        78 -99 125 -123 -111 -71 17
        -115 78 40 -88 81 -40 78
        -99 126 -54 1 0 0 24
        -55 88 42 0 1 0 -8
        119 55 2 0 0 1 -92
        -40 -21 -89 -79 59 -44 9
    ],) == true
    @test has_identity(
        [
            78 -99 125 -123 -111 -71 17
            -115 78 40 -88 81 -40 78
            -99 126 -54 1 0 0 24
            -55 88 42 0 1 0 -8
            119 55 2 0 0 1 -92
            -40 -21 -89 -79 59 -44 9
        ],
        4,
    ) == false
    @test has_identity_on_left([1 0 0 2 3 0; 0 1 0 1 2 2; 0 0 1 4 3 0]) == true
    @test has_identity_on_left([-96 -66 20 1 0 0; -65 59 -82 0 1 0; -16 87 -113 0 0 1]) ==
        false

    @test hamming_distance("ABC", "DBC") == 1
    @test hamming_distance("ABC", "DEF") == 3

    @test aresamearrays(
        hamming_ball([[1, 0, 1], [0, 1, 1], [1, 0, 0]], [1, 0, 0], 2),
        deepsym([[1, 0, 1], [1, 0, 0]]),
    )
    @test aresamearrays(
        hamming_ball(list_span([2, 2, 2], [1, 2, 0], [0, 2, 1], 3), [1, 0, 0], 3),
        deepsym([
            [0, 0, 0],
            [0, 2, 1],
            [0, 1, 2],
            [1, 2, 0],
            [1, 1, 1],
            [1, 0, 2],
            [2, 1, 0],
            [2, 0, 1],
            [2, 2, 2],
        ]),
    )

    @test rate(3, 5, 4) ≈ 0.3662433802

    @test t_error_correcting([[0, 0, 0, 0], [2, 2, 2, 2]], 1) == true
    @test t_error_correcting([[1, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 1]], 3) == false
    @test t_error_detecting([[1, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 1]], 3) == false

    @test find_error_detection_max(
        [[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 1, 0], [1, 1, 0, 1]], 2
    ) == 1
    @test find_error_correction_max(
        [[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 1, 0], [1, 1, 0, 1]], 2
    ) == 0
    @test find_error_correction_max(list_span([1, 0, 1, 0], [0, 1, 1, 1], 2), 2) == 0
    @test find_error_detection_max(list_span([1, 0, 1, 0], [0, 1, 1, 1], 2), 2) == 1

    @test Polynomial([1, 2, 3, 4, 5, 6, 7, 8, 9], 3) == Polynomial([1, 2, 0, 1, 2, 0, 1, 2])

    @test isirreducible(Polynomial([1, 1, 0, 0, 1]), 2) == true
    @test isirreducible(Polynomial([1, 1, 1, 0, 1, 1]), 2) == true
    @test isirreducible(Polynomial([4, 4, 1]), 2) == false
    @test isirreducible(Polynomial([1, 0, 1]), 2) == false
    @test isirreducible(Polynomial([-2, 0, 1]), 2) == false
    @test isirreducible(Polynomial([1, 1]), 2) == false

    p = Polynomial([1, 1, 2, 0, 1, 2, 1])
    q = Polynomial([2, 1, 1])
    a = Polynomial([0, 1, 0, 1, 0, 0, 1])
    b = Polynomial([1, 0, 1])

    @test mod(rem(p, q), 3) == Polynomial([1, 1])
    @test mod(rem(a, b), 2) == Polynomial([1])

    @test isperfectpower(36) == true
    @test isprimepower(36) == false
    @test isperfectpower(9) == true
    @test isprimepower(9) == true
    @test isperfectpower(5) == false
    @test isprimepower(2) == true
    @test isprimepower(18) == false

    @test rref([1 0 1 1 1; 1 1 1 0 1; 0 1 1 1 1], 2) == [1 0 0 1 0; 0 1 0 1 0; 0 0 1 0 1]
    @test rref([1 0 1 0 1 0; 0 1 0 0 1 0; 1 1 1 1 1 1], 2) ==
        [1 0 1 0 1 0; 0 1 0 0 1 0; 0 0 0 1 1 1]
    @test rref([1 1 0 2 3 1; 2 0 1 3 4 1; 1 2 2 1 4 3], 5; colswap = false) ==
        [1 0 3 0 2 2; 0 1 2 0 1 1; 0 0 0 1 0 4]
    @test rref([1 1 0 2 3 1; 2 0 1 3 4 1; 1 2 2 1 4 3], 5; colswap = true) ==
        [1 0 0 3 2 2; 0 1 0 2 1 1; 0 0 1 0 0 4]
    @test rref([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) ==
        [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1]
    @test rref([0 0 0 0 0; 1 0 1 0 1; 0 1 0 1 1; 1 1 1 1 0], 2) ==
        [1 0 1 0 1; 0 1 0 1 1; 0 0 0 0 0; 0 0 0 0 0]
    @test rref([1 1 1 0; 1 1 0 1; 0 0 1 1], 2) == [1 1 0 1; 0 0 1 1; 0 0 0 0]

    @test multiplication_table(2, 3) == Polynomial[
        Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0])
        Polynomial([0]) Polynomial([1]) Polynomial([2]) Polynomial([0, 1]) Polynomial([1, 1]) Polynomial([2, 1]) Polynomial([0, 2]) Polynomial([1, 2]) Polynomial([2, 2])
        Polynomial([0]) Polynomial([2]) Polynomial([1]) Polynomial([0, 2]) Polynomial([2, 2]) Polynomial([1, 2]) Polynomial([0, 1]) Polynomial([2, 1]) Polynomial([1, 1])
        Polynomial([0]) Polynomial([0, 1]) Polynomial([0, 2]) Polynomial([0, 0, 1]) Polynomial([0, 1, 1]) Polynomial([0, 2, 1]) Polynomial([0, 0, 2]) Polynomial([0, 1, 2]) Polynomial([0, 2, 2])
        Polynomial([0]) Polynomial([1, 1]) Polynomial([2, 2]) Polynomial([0, 1, 1]) Polynomial([1, 2, 1]) Polynomial([2, 0, 1]) Polynomial([0, 2, 2]) Polynomial([1, 0, 2]) Polynomial([2, 1, 2])
        Polynomial([0]) Polynomial([2, 1]) Polynomial([1, 2]) Polynomial([0, 2, 1]) Polynomial([2, 0, 1]) Polynomial([1, 1, 1]) Polynomial([0, 1, 2]) Polynomial([2, 2, 2]) Polynomial([1, 0, 2])
        Polynomial([0]) Polynomial([0, 2]) Polynomial([0, 1]) Polynomial([0, 0, 2]) Polynomial([0, 2, 2]) Polynomial([0, 1, 2]) Polynomial([0, 0, 1]) Polynomial([0, 2, 1]) Polynomial([0, 1, 1])
        Polynomial([0]) Polynomial([1, 2]) Polynomial([2, 1]) Polynomial([0, 1, 2]) Polynomial([1, 0, 2]) Polynomial([2, 2, 2]) Polynomial([0, 2, 1]) Polynomial([1, 1, 1]) Polynomial([2, 0, 1])
        Polynomial([0]) Polynomial([2, 2]) Polynomial([1, 1]) Polynomial([0, 2, 2]) Polynomial([2, 1, 2]) Polynomial([1, 0, 2]) Polynomial([0, 1, 1]) Polynomial([2, 0, 1]) Polynomial([1, 2, 1])
    ]

    @test aresamearrays(
        list_span([2, 1, 1], [1, 1, 1], 3),
        [
            [0, 0, 0],
            [1, 1, 1],
            [2, 2, 2],
            [2, 1, 1],
            [0, 2, 2],
            [1, 0, 0],
            [1, 2, 2],
            [2, 0, 0],
            [0, 1, 1],
        ],
    )

    @test islinear([[0, 0, 0], [1, 1, 1], [1, 0, 1], [1, 1, 0]], 2) == false
    @test islinear([[0, 0, 0], [1, 1, 1], [1, 0, 1], [0, 1, 0]], 2) == true

    @test code_distance([
        [0, 0, 0, 0, 0], [1, 0, 1, 0, 1], [0, 1, 0, 1, 0], [1, 1, 1, 1, 1]
    ]) == 2
    @test code_distance([
        [0, 0, 0, 0, 0],
        [1, 1, 1, 0, 0],
        [0, 0, 0, 1, 1],
        [1, 1, 1, 1, 1],
        [1, 0, 0, 1, 1],
        [0, 1, 1, 0, 0],
    ]) == 1

    w1 = Word((1, 5, 2, 6))
    w2 = Word([4, 2, 6, 7, 3, 10, 10])
    w3 = Word(1, 3, 6, 2, 4)
    @test all(isword.([w1, w2, w3])) == true
    @test all(isabstractword.([w1, w2, w3])) == true
    @test isabstractword(1, 2, 3) == true
    @test all(isabstractword.([(1, 2, 3), [1, 2, 3]])) == true
    @test any(isword.([(1, 2, 3), [4, 5, 6]])) == false
    @test isword(1, 2, 3) == false
    @test aresamearrays(Alphabet("123"), deepsym([1, 2, 3]))
    @test aresamearrays(Alphabet([1, 2, 3]), deepsym([1, 2, 3]))
    @test aresamearrays(Alphabet(["1", "2", "3"]), deepsym([1, 2, 3]))
    # TODO: write test for CodeUniverse struct
    @test aresamearrays(
        [i for i in CodeUniverseIterator(["a", "b", "c"], 3)],
        get_all_words(["a", "b", "c"], 3),
    ) # implicitly tests CodeUniverseIterator
    @test aresamearrays(
        collect(CodeUniverseIterator(["a", "b", "c"], 4)),
        Tuple[
            (:a, :a, :a, :a),
            (:b, :a, :a, :a),
            (:c, :a, :a, :a),
            (:a, :b, :a, :a),
            (:b, :b, :a, :a),
            (:c, :b, :a, :a),
            (:a, :c, :a, :a),
            (:b, :c, :a, :a),
            (:c, :c, :a, :a),
            (:a, :a, :b, :a),
            (:b, :a, :b, :a),
            (:c, :a, :b, :a),
            (:a, :b, :b, :a),
            (:b, :b, :b, :a),
            (:c, :b, :b, :a),
            (:a, :c, :b, :a),
            (:b, :c, :b, :a),
            (:c, :c, :b, :a),
            (:a, :a, :c, :a),
            (:b, :a, :c, :a),
            (:c, :a, :c, :a),
            (:a, :b, :c, :a),
            (:b, :b, :c, :a),
            (:c, :b, :c, :a),
            (:a, :c, :c, :a),
            (:b, :c, :c, :a),
            (:c, :c, :c, :a),
            (:a, :a, :a, :b),
            (:b, :a, :a, :b),
            (:c, :a, :a, :b),
            (:a, :b, :a, :b),
            (:b, :b, :a, :b),
            (:c, :b, :a, :b),
            (:a, :c, :a, :b),
            (:b, :c, :a, :b),
            (:c, :c, :a, :b),
            (:a, :a, :b, :b),
            (:b, :a, :b, :b),
            (:c, :a, :b, :b),
            (:a, :b, :b, :b),
            (:b, :b, :b, :b),
            (:c, :b, :b, :b),
            (:a, :c, :b, :b),
            (:b, :c, :b, :b),
            (:c, :c, :b, :b),
            (:a, :a, :c, :b),
            (:b, :a, :c, :b),
            (:c, :a, :c, :b),
            (:a, :b, :c, :b),
            (:b, :b, :c, :b),
            (:c, :b, :c, :b),
            (:a, :c, :c, :b),
            (:b, :c, :c, :b),
            (:c, :c, :c, :b),
            (:a, :a, :a, :c),
            (:b, :a, :a, :c),
            (:c, :a, :a, :c),
            (:a, :b, :a, :c),
            (:b, :b, :a, :c),
            (:c, :b, :a, :c),
            (:a, :c, :a, :c),
            (:b, :c, :a, :c),
            (:c, :c, :a, :c),
            (:a, :a, :b, :c),
            (:b, :a, :b, :c),
            (:c, :a, :b, :c),
            (:a, :b, :b, :c),
            (:b, :b, :b, :c),
            (:c, :b, :b, :c),
            (:a, :c, :b, :c),
            (:b, :c, :b, :c),
            (:c, :c, :b, :c),
            (:a, :a, :c, :c),
            (:b, :a, :c, :c),
            (:c, :a, :c, :c),
            (:a, :b, :c, :c),
            (:b, :b, :c, :c),
            (:c, :b, :c, :c),
            (:a, :c, :c, :c),
            (:b, :c, :c, :c),
            (:c, :c, :c, :c),
        ],
    )

    @test sphere_covering_bound(5, 7, 3) == 215
    @test sphere_packing_bound(5, 7, 3) == 2693

    @test aresamearrays(
        construct_ham_matrix(3, 2), [0 0 0 1 1 1 1; 0 1 1 0 0 1 1; 1 0 1 0 1 0 1]
    )
    @test aresamearrays(
        construct_ham_matrix(3, 3),
        [0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 1 1 1 2 2 2 0 0 0 1 1; 1 2 0 1 2 0 1 2 0 1 2 0 1],
    )

    @test isperfect(11, 6, 5, 3) == true
    @test isperfect(23, 12, 7, 2) == true
    @test isperfect(23, 12, 7, 3) == false
    @test isperfect(11, 6, 5, 4) == false
    @test isgolayperfect(11, 6, 5, 3) == true
    @test isgolayperfect(23, 12, 7, 2) == true
    @test isgolayperfect(23, 12, 7, 3) == false
    @test isgolayperfect(11, 6, 5, 4) == false

    @test length(get_codewords(5, 5, 3)) ∈ [74:74...] # implicitly tests Word struct, and `get_codewords_random`
    @test length(get_codewords(4, 7, 3; m = 1)) ∈ [256:308...]
    @test length(get_codewords_greedy(5, 5, 3)) == 74
    randq, randn = rand(1:8, 2)
    @test length(get_all_words(randq, randn)) == big(randq)^randn
    @test aresamearrays(
        get_codewords([1 0 1 0; 0 1 1 1], 2),
        [[0, 0, 0, 0], [1, 0, 1, 0], [0, 1, 1, 1], [1, 1, 0, 1]],
    )
    @test aresamearrays(
        get_codewords([1 0 0 1 1 0; 0 1 0 1 0 1; 0 0 1 0 1 1], 2),
        [
            [0, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 1, 0],
            [0, 1, 0, 1, 0, 1],
            [1, 1, 0, 0, 1, 1],
            [0, 0, 1, 0, 1, 1],
            [1, 0, 1, 1, 0, 1],
            [0, 1, 1, 1, 1, 0],
            [1, 1, 1, 0, 0, 0],
        ],
    )

    # note that the following tests are not to be checked using `aresamearrays` because their order matters
    @test syndrome(
        [0, 2, 1, 2, 0, 1, 0],
        transpose(parity_check(
            [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3
        )),
        3,
    ) == [0 0 0]
    @test parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3) ==
        [1 1 2 1 1 0 0; 1 0 0 1 0 1 0; 1 2 1 2 0 0 1]
    @test normal_form([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) ==
        [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1]
    @test equivalent_code(
        [1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3
    ) == [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1]
    @test parity_check(
        normal_form([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3), 3
    ) == [1 1 2 1 1 0 0; 1 0 0 1 0 1 0; 1 2 1 2 0 0 1]
    @test isincode(
        [0, 2, 1, 2, 0, 1, 0],
        transpose(parity_check(
            [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3
        )),
        3,
    ) == true
    @test isincode(
        [1, 0, 2, 2, 1, 2, 1],
        transpose(parity_check(
            [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3
        )),
        3,
    ) == false
end # end runtests

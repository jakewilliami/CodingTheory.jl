using CodingTheory, BenchmarkTools, Polynomials
import CodingTheory.deepsym, CodingTheory.displaymatrix

function BM()
    hamming_distance("ABC", "DBC") == 1
    hamming_distance("ABC", "DEF") == 3

    hamming_ball([[1, 0, 1], [0, 1, 1], [1, 0, 0]], [1, 0, 0], 2) == deepsym([[1, 0, 1], [1, 0, 0]])
    hamming_ball(list_span([2, 2, 2], [1, 2, 0], [0, 2, 1], 3), [1, 0, 0], 3) == deepsym([[0, 0, 0], [0, 2, 1], [0, 1, 2], [1, 2, 0], [1, 1, 1], [1, 0, 2], [2, 1, 0], [2, 0, 1], [2, 2, 2]])

    rate(3, 5, 4) ≈ 0.3662433802

    t_error_correcting([[0, 0, 0, 0], [2, 2, 2, 2]], 1) == true
    t_error_correcting([[1, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 1]], 3) == false
    t_error_detecting([[1, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 1]], 3) == false

    find_error_detection_max([[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 1, 0], [1, 1, 0, 1]], 2) == 1
    find_error_correction_max([[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 1, 0], [1, 1, 0, 1]], 2) == 0
    find_error_correction_max(list_span([1, 0, 1, 0], [0, 1, 1, 1], 2), 2) == 0
    find_error_detection_max(list_span([1, 0, 1, 0], [0, 1, 1, 1], 2), 2) == 1

    Polynomial([1, 2, 3, 4, 5, 6, 7, 8, 9], 3) == Polynomial([1, 2, 0, 1, 2, 0, 1, 2])

    isirreducible(Polynomial([1, 1, 0, 0, 1]), 2) == true
    isirreducible(Polynomial([1, 1, 1, 0, 1, 1]), 2) == true
    isirreducible(Polynomial([4, 4, 1]), 2) == false
    isirreducible(Polynomial([1, 0, 1]), 2) == false
    isirreducible(Polynomial([-2, 0, 1]), 2) == false
    isirreducible(Polynomial([1, 1]), 2) == false

    p = Polynomial([1, 1, 2, 0, 1, 2, 1])
    q = Polynomial([2, 1, 1])
    a = Polynomial([0, 1, 0, 1, 0, 0, 1])
    b = Polynomial([1, 0, 1])

    mod(rem(p, q), 3) == Polynomial([1, 1])
    mod(rem(a, b), 2) == Polynomial([1])

    rref([1 0 1 1 1; 1 1 1 0 1; 0 1 1 1 1], 2) == [1 0 0 1 0; 0 1 0 1 0; 0 0 1 0 1]
    rref([1 0 1 0 1 0; 0 1 0 0 1 0; 1 1 1 1 1 1], 2) == [1 0 1 0 1 0; 0 1 0 0 1 0; 0 0 0 1 1 1]
    rref([1 1 0 2 3 1; 2 0 1 3 4 1; 1 2 2 1 4 3], 5, colswap=false) == [1 0 3 0 2 2; 0 1 2 0 1 1; 0 0 0 1 0 4]
    rref([1 1 0 2 3 1; 2 0 1 3 4 1; 1 2 2 1 4 3], 5, colswap=true) == [1 0 0 3 2 2; 0 1 0 2 1 1; 0 0 1 0 0 4]
    rref([1 2 0 1 2 1 2;  2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) == [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1]
    rref([0 0 0 0 0; 1 0 1 0 1; 0 1 0 1 1; 1 1 1 1 0], 2) == [1 0 1 0 1; 0 1 0 1 1; 0 0 0 0 0; 0 0 0 0 0]
    rref([1 1 1 0; 1 1 0 1; 0 0 1 1], 2) == [1 1 0 1; 0 0 1 1; 0 0 0 0]

    multiplication_table(2, 3) == Polynomial[Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]) Polynomial([0]); Polynomial([0]) Polynomial([1]) Polynomial([2]) Polynomial([0, 1]) Polynomial([1, 1]) Polynomial([2, 1]) Polynomial([0, 2]) Polynomial([1, 2]) Polynomial([2, 2]); Polynomial([0]) Polynomial([2]) Polynomial([1]) Polynomial([0, 2]) Polynomial([2, 2]) Polynomial([1, 2]) Polynomial([0, 1]) Polynomial([2, 1]) Polynomial([1, 1]); Polynomial([0]) Polynomial([0, 1]) Polynomial([0, 2]) Polynomial([0, 0, 1]) Polynomial([0, 1, 1]) Polynomial([0, 2, 1]) Polynomial([0, 0, 2]) Polynomial([0, 1, 2]) Polynomial([0, 2, 2]); Polynomial([0]) Polynomial([1, 1]) Polynomial([2, 2]) Polynomial([0, 1, 1]) Polynomial([1, 2, 1]) Polynomial([2, 0, 1]) Polynomial([0, 2, 2]) Polynomial([1, 0, 2]) Polynomial([2, 1, 2]); Polynomial([0]) Polynomial([2, 1]) Polynomial([1, 2]) Polynomial([0, 2, 1]) Polynomial([2, 0, 1]) Polynomial([1, 1, 1]) Polynomial([0, 1, 2]) Polynomial([2, 2, 2]) Polynomial([1, 0, 2]); Polynomial([0]) Polynomial([0, 2]) Polynomial([0, 1]) Polynomial([0, 0, 2]) Polynomial([0, 2, 2]) Polynomial([0, 1, 2]) Polynomial([0, 0, 1]) Polynomial([0, 2, 1]) Polynomial([0, 1, 1]); Polynomial([0]) Polynomial([1, 2]) Polynomial([2, 1]) Polynomial([0, 1, 2]) Polynomial([1, 0, 2]) Polynomial([2, 2, 2]) Polynomial([0, 2, 1]) Polynomial([1, 1, 1]) Polynomial([2, 0, 1]); Polynomial([0]) Polynomial([2, 2]) Polynomial([1, 1]) Polynomial([0, 2, 2]) Polynomial([2, 1, 2]) Polynomial([1, 0, 2]) Polynomial([0, 1, 1]) Polynomial([2, 0, 1]) Polynomial([1, 2, 1])]

    list_span([2, 1, 1], [1, 1, 1], 3) == [[0, 0, 0], [1, 1, 1], [2, 2, 2], [2, 1, 1], [0, 2, 2], [1, 0, 0], [1, 2, 2], [2, 0, 0], [0, 1, 1]]

    islinear([[0,0,0],[1,1,1],[1,0,1],[1,1,0]], 2) == false
    islinear([[0,0,0],[1,1,1],[1,0,1],[0,1,0]], 2) == true

    code_distance([[0,0,0,0,0],[1,0,1,0,1],[0,1,0,1,0],[1,1,1,1,1]]) == 2
    code_distance([[0,0,0,0,0],[1,1,1,0,0],[0,0,0,1,1],[1,1,1,1,1],[1,0,0,1,1],[0,1,1,0,0]]) == 1

    Alphabet("123") == deepsym([1, 2, 3])
    Alphabet([1, 2, 3]) == deepsym([1, 2, 3])
    Alphabet(["1", "2", "3"]) == deepsym([1, 2, 3])
    # TODO: write test for CodeUniverse struct
    [i for i in CodeUniverseIterator(["a", "b", "c"], 3)] == get_all_words(["a", "b", "c"], 3) # implicitly tests CodeUniverseIterator
    collect(CodeUniverseIterator(["a", "b", "c"], 4)) == Tuple[(:a, :a, :a, :a), (:b, :a, :a, :a), (:c, :a, :a, :a), (:a, :b, :a, :a), (:b, :b, :a, :a), (:c, :b, :a, :a), (:a, :c, :a, :a), (:b, :c, :a, :a), (:c, :c, :a, :a), (:a, :a, :b, :a), (:b, :a, :b, :a), (:c, :a, :b, :a), (:a, :b, :b, :a), (:b, :b, :b, :a), (:c, :b, :b, :a), (:a, :c, :b, :a), (:b, :c, :b, :a), (:c, :c, :b, :a), (:a, :a, :c, :a), (:b, :a, :c, :a), (:c, :a, :c, :a), (:a, :b, :c, :a), (:b, :b, :c, :a), (:c, :b, :c, :a), (:a, :c, :c, :a), (:b, :c, :c, :a), (:c, :c, :c, :a), (:a, :a, :a, :b), (:b, :a, :a, :b), (:c, :a, :a, :b), (:a, :b, :a, :b), (:b, :b, :a, :b), (:c, :b, :a, :b), (:a, :c, :a, :b), (:b, :c, :a, :b), (:c, :c, :a, :b), (:a, :a, :b, :b), (:b, :a, :b, :b), (:c, :a, :b, :b), (:a, :b, :b, :b), (:b, :b, :b, :b), (:c, :b, :b, :b), (:a, :c, :b, :b), (:b, :c, :b, :b), (:c, :c, :b, :b), (:a, :a, :c, :b), (:b, :a, :c, :b), (:c, :a, :c, :b), (:a, :b, :c, :b), (:b, :b, :c, :b), (:c, :b, :c, :b), (:a, :c, :c, :b), (:b, :c, :c, :b), (:c, :c, :c, :b), (:a, :a, :a, :c), (:b, :a, :a, :c), (:c, :a, :a, :c), (:a, :b, :a, :c), (:b, :b, :a, :c), (:c, :b, :a, :c), (:a, :c, :a, :c), (:b, :c, :a, :c), (:c, :c, :a, :c), (:a, :a, :b, :c), (:b, :a, :b, :c), (:c, :a, :b, :c), (:a, :b, :b, :c), (:b, :b, :b, :c), (:c, :b, :b, :c), (:a, :c, :b, :c), (:b, :c, :b, :c), (:c, :c, :b, :c), (:a, :a, :c, :c), (:b, :a, :c, :c), (:c, :a, :c, :c), (:a, :b, :c, :c), (:b, :b, :c, :c), (:c, :b, :c, :c), (:a, :c, :c, :c), (:b, :c, :c, :c), (:c, :c, :c, :c)]

    sphere_covering_bound(5,7,3) == 215
    sphere_packing_bound(5,7,3) == 2693

    construct_ham_matrix(3,2) == [0 0 0 1 1 1 1; 0 1 1 0 0 1 1; 1 0 1 0 1 0 1]
    construct_ham_matrix(3,3) == [0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 1 1 1 2 2 2 0 0 0 1 1; 1 2 0 1 2 0 1 2 0 1 2 0 1]

    isperfect(11, 6, 5, 3) == true
    isperfect(23, 12, 7, 2) == true
    isperfect(23, 12, 7, 3) == false
    isperfect(11, 6, 5, 4) == false
    isgolayperfect(11, 6, 5, 3) == true
    isgolayperfect(23, 12, 7, 2) == true
    isgolayperfect(23, 12, 7, 3) == false
    isgolayperfect(11, 6, 5, 4) == false

    length(get_codewords(5, 5, 3)) ∈ [74:74...]
    length(get_codewords(4, 7, 3; m = 1)) ∈ [256:308...]
    length(get_codewords_greedy(5, 5, 3)) == 74
    randq, randn = rand(1:8, 2)
    length(get_all_words(randq, randn)) == big(randq)^randn
    get_codewords([1 0 1 0; 0 1 1 1], 2) == [[0, 0, 0, 0], [1, 0, 1, 0], [0, 1, 1, 1], [1, 1, 0, 1]]
    get_codewords([1 0 0 1 1 0; 0 1 0 1 0 1; 0 0 1 0 1 1], 2) == [[0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 1, 0], [0, 1, 0, 1, 0, 1], [1, 1, 0, 0, 1, 1], [0, 0, 1, 0, 1, 1], [1, 0, 1, 1, 0, 1], [0, 1, 1, 1, 1, 0], [1, 1, 1, 0, 0, 0]]

    syndrome([0, 2, 1, 2, 0, 1, 0], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3) == [0 0 0]
    parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3) == [1 1 2 1 1 0 0; 1 0 0 1 0 1 0; 1 2 1 2 0 0 1]
    normal_form([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) == [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1]
    equivalent_code([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) == [1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1]
    parity_check(normal_form([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3), 3) == [1 1 2 1 1 0 0; 1 0 0 1 0 1 0; 1 2 1 2 0 0 1]
    isincode([0, 2, 1, 2, 0, 1, 0], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3) == true
    isincode([1, 0, 2, 2, 1, 2, 1], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3) == false
end # end runtests

@btime BM()

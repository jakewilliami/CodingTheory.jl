<h1 align="center">
    CodingTheory.jl
</h1>

[![Code Style: Blue][code-style-img]][code-style-url] [![Build Status](https://travis-ci.com/jakewilliami/CodingTheory.jl.svg?branch=master)](https://travis-ci.com/jakewilliami/CodingTheory.jl) ![Project Status](https://img.shields.io/badge/status-maturing-green)

## Description
This is a minimal package for a pure Julia implementation of tools used in [Coding Theory](https://en.wikipedia.org/wiki/Coding_theory).  The is the science of accurately transmitting information through a noisy channel.

## Terminology
We assume that Alice and Bob communicate by sending sequences of symbols from a finite set *&Sigma;*, which we call the **alphabet**.  We always use *q* to stand for the size of the set of symbols, |*&Sigma;*|.  A **word** is a sequence of symbols from the alphabet *&Sigma;*.  If *w<sub>1</sub>w<sub>2</sub>...w<sub>n</sub>* is such a word, then *n* is its **length**.  We use *&Sigma;<sup>n</sup>* to denote the set of words with length *n* using symbols in *&Sigma;*.  In general, the number of words in *&Sigma;<sup>n</sup>* is
<p align="center">
    |&Sigma;<sup>n</sup>| = q<sup>n</sup>.
</p>

**Block codes** are codes in which Alice transmits words of a preditermined and fixed length.  A **code** is a subset *C &SubsetEqual; &Sigma;<sup>n</sup>*.  The words in *C* are called **code words**.  We say that *n* is the **block length**.  We use *M* to stand for *|C|*, the number of code words.  Alice has a set *&Mellintrf;*, some of which she wants to send to Bob, so she has the bijective encoding function
<p align="center">
    E : &Mellintrf; &longrightarrow; C.
</p>  

Similarly, Bob has a decoding function
<p align="center">
    D : &Sigma;<sup>n</sup> &longrightarrow; C &cup; {?},
</p>

Where Bob uses the *?* symbol when he cannot confidently decode.  So if Alice wishes to communicate a message, she transmits a code word *w = E(M)*.  *w* may be corrupted to *w' &ne; w*.  Then Bob can decode *w'* as *E<sup>-1</sup>(D(w))*.  If Bob is not certain how to decode, then *D(w')* may be '*?*', which means that Bob can tell an error has occurred but is not certain what that error is.

If *&Mellintrf; &SubsetEqual; &Sigma;<sup>k</sup>* is the set of messages, then *k* is the **message length**.


## Examples

```julia
julia> using CodingTheory

julia> hamming_distance("ABC", "BBC") # computes the hamming distance
1

julia> hamming_ball([[1, 0, 1], [0, 1, 1], [1, 0, 0]], [1, 0, 0], 2) # given a list of words, a word, and a distance e (respectively), calculate all the words in the alphabet within distance e of that word.  Converts to symbols in order to keep unique lengths
2-element Array{Any,1}:
 [Symbol("1"), Symbol("0"), Symbol("1")]
 [Symbol("1"), Symbol("0"), Symbol("0")]

julia> t_error_detecting([[1, 0, 1], [0, 1, 1], [1, 0, 0]], 3)
false

julia> find_error_detection_max([[0, 0, 0, 0], [0, 1, 1, 1], [1, 0, 1, 0], [1, 1, 0, 1]], 2)
1

julia> isirreducible(Polynomial([1, 1, 0, 0, 1]), 2) # is 1 + x + x^4 mod 2 irreducible?
true

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

julia> islinear([[0,0,0],[1,1,1],[1,0,1],[1,1,0]], 2) # checks whether a vector of vectors is linear/a subspace (modulo 2)
false

julia> code_distance([[0,0,0,0,0],[1,0,1,0,1],[0,1,0,1,0],[1,1,1,1,1]]) # gets the minimum distance between two vectors in an array of vectors
2

julia> rate(3, 5, 4) # the rate of the code which has 3 symbols, 5 words in the code, and word length of 4 (e.g., Σ = {A, B, C}, C = {ABBA,CABA,BBBB,CAAB,ACBB})
0.3662433801794817

julia> sphere_covering_bound(5,7,3)
215

julia> sphere_packing_bound(5,7,3)
2693

julia> construct_ham_matrix(3,2)
3×7 Array{Int64,2}:
 0  0  0  1  1  1  1
 0  1  1  0  0  1  1
 1  0  1  0  1  0  1

julia> isperfect(11, 6, 5, 3)
true

julia> isgolayperfect(11, 6, 5, 3)
true

julia> ishammingperfect(11, 6, 5, 3)
false

julia> rref([1 1 0 2 3 1; 2 0 1 3 4 1; 1 2 2 1 4 3], 5, colswap=true) # gauss-jordan elimitation modulo 5 with column swapping
3×6 Array{Int64,2}:
 1  0  0  3  2  2
 0  1  0  2  1  1
 0  0  1  0  0  4

julia> get_codewords(["a", "b", "c"], 3, 2) # get codewords of block length 3 with distance 2.  Once again, are symbols for uniqueness
9-element Array{Tuple,1}:
 (:c, :b, :b)
 (:b, :c, :b)
 (:b, :b, :c)
 (:c, :c, :a)
 (:c, :a, :c)
 (:a, :c, :c)
 (:a, :b, :a)
 (:b, :a, :a)
 (:a, :a, :b)

julia> get_all_words(2, 2) # all words of block length 2 using 2 unique symbols
4-element Array{Tuple,1}:
(Symbol("##258"), Symbol("##258"))
(Symbol("##259"), Symbol("##258"))
(Symbol("##258"), Symbol("##259"))
(Symbol("##259"), Symbol("##259"))

julia> syndrome([0, 2, 1, 2, 0, 1, 0], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3)
1×3 Array{Int64,2}:
 0  0  0

julia> normal_form([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) # computes rref colswap = false
4×7 Array{Int64,2}:
1  0  0  0  2  2  2
0  1  0  0  2  0  1
0  0  1  0  1  0  2
0  0  0  1  2  2  1

julia> equivalent_code([1 2 0 1 2 1 2; 2 2 2 0 1 1 1; 1 0 1 1 2 1 2; 0 1 0 1 1 2 2], 3) # computes rref colswap = true
4×7 Array{Int64,2}:
1  0  0  0  2  2  2
0  1  0  0  2  0  1
0  0  1  0  1  0  2
0  0  0  1  2  2  1

julia> isincode([0, 2, 1, 2, 0, 1, 0], transpose(parity_check([1 0 0 0 2 2 2; 0 1 0 0 2 0 1; 0 0 1 0 1 0 2; 0 0 0 1 2 2 1], 3)), 3) # tests if the syndrome is equal to the zero vector, and is thus in the code
true
```

## A note on the number of codewords in a code

We have some algorithms brute-force searching for the codewords in a [q, n, d]-code.  These algorithms are brute-force as they do not assume that q is a prime power.  Therefore, they go through all possible codewords of a [q, n]-code, and narrow down the code based on d.  There algorithms are namely `get_codewords_greedy` and `get_codewords_random`, both of which using `get_all_codewords`.  The `get_codewords` function iterates through possibilities of `get_codeword_random` and chooses the maximum of those iterations or the `get_codeword_greedy` length.  Despite the name, `get_codewords` is only a **probably** candidate.  Increate the keyword argument `m` to decrease the likelihood that there is a code with more codewords while maintaining the bound of the distance.  Furthermore, there is a `get_codewords` method that lists all linear combinations of rows of a generator matrix.

[code-style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[code-style-url]: https://github.com/invenia/BlueStyle

<h1 align="center">
    CodingTheory.jl
</h1>

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jakewilliami.github.io/CodingTheory.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jakewilliami.github.io/CodingTheory.jl/dev)
[![CI](https://github.com/invenia/PkgTemplates.jl/workflows/CI/badge.svg)](https://github.com/jakewilliami/CodingTheory.jl/actions?query=workflow%3ACI)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
![Project Status](https://img.shields.io/badge/status-maturing-green)

This is a minimal package for a pure Julia implementation of tools used in [Coding Theory](https://en.wikipedia.org/wiki/Coding_theory).  This is the science of accurately transmitting information through a noisy channel.

## Background
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

## A note on the number of codewords in a code

We have some algorithms brute-force searching for the codewords in a [q, n, d]-code.  These algorithms are brute-force as they do not assume that q is a prime power.  Therefore, they go through all possible codewords of a [q, n]-code, and narrow down the code based on d.  There algorithms are namely `get_codewords_greedy` and `get_codewords_random`, both of which using `get_all_codewords`.  The `get_codewords` function iterates through possibilities of `get_codeword_random` and chooses the maximum of those iterations or the `get_codeword_greedy` length.  Despite the name, `get_codewords` is only a **probably** candidate.  Increate the keyword argument `m` to decrease the likelihood that there is a code with more codewords while maintaining the bound of the distance.  Furthermore, there is a `get_codewords` method that lists all linear combinations of rows of a generator matrix.


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


[code-style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[code-style-url]: https://github.com/invenia/BlueStyle

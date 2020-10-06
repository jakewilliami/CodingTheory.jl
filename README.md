<h1 align="center">
    CodingTheory.jl
</h1>


## Description
This is a minimal package for a pure Julia implementation of tools used in [Coding Theory](https://en.wikipedia.org/wiki/Coding_theory).  The is the science of accurately transmitting information through a noisy channel.

## Terminology
We assume that Alice and Bob communicate by sending sequences of symbols from a finite set *&Sigma;*, which we call the **alphabet**.  We always use *q* to stand for the size of the set of symbols, |*&Sigma;*|.  A **word** is a sequence of symbols from the alphabet *&Sigma;*.  If *w<sub>1</sub>w<sub>2</sub>...w<sub>n</sub>* is such a word, then *n* is its **length**.  We use *&Sigma;<sup>n</sup>* to denote the set of words with length *n* using symbols in *&Sigma;*.  In general, the number of words in *&Sigma;<sup>n</sup>* is
<p align="center">
    |&Sigma;<sup>n</sup>|=q<sup>n</sup>.
</p>

**Block codes** are codes in which Alice transmits words of a preditermined and fixed length.  A **code** is a subset *C&SubsetEqual;&Sigma;<sup>n</sup>*.

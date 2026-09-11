# Mutually unbiased bases in dimension six (`mubs`)

Two orthonormal bases {e_i}, {f_j} of C^d are **mutually unbiased** (MU)
if |<e_i|f_j>|^2 = 1/d for all i, j. A set of MUBs is a set of bases
that are pairwise mutually unbiased. In every prime-power dimension d,
exactly d + 1 MUBs exist. Dimension six is the first dimension where the
answer is unknown; **Zauner's conjecture** (1999) states that no more
than three MUBs exist in C^6, far below the upper bound of seven. This
is one of the best-known open problems in quantum information theory.

Any set of MUBs can be written as {I, H1/sqrt(6), H2/sqrt(6), ...} where
the Hi are complex Hadamard matrices (unimodular entries, orthogonal
columns). The conjecture is therefore equivalent to: no three complex
Hadamard matrices of order six are pairwise mutually unbiased and
unbiased to I (no "MU quadruple" of bases). Progress typically takes the
form of excluding specific strata of the complex Hadamard landscape from
participating in a quadruple.

**Butson strata.** A Butson Hadamard matrix BH(d, N) is a complex
Hadamard matrix whose entries are all N-th roots of unity. These
discrete strata admit finite, exact-arithmetic exhaustion, and they are
where all known complete MUB sets in prime-power dimensions live, which
makes them natural targets.

## Claims in this problem

| claim | statement (short) | status |
|---|---|---|
| [claim-164](claims/claim-164-butson-n24/) | No MU quadruple in C^6 from Butson Hadamard matrices of root order 24 (hence any root order dividing 24, and, via a corollary, any single root order N <= 12). | internally certified, awaiting external verification |

## Standard references

- I. Bengtsson, W. Bruzda, A. Ericsson, J.-A. Larsson, W. Tadej,
  K. Zyczkowski, "Mutually unbiased bases and Hadamard matrices of
  order six", J. Math. Phys. 48, 052106 (2007), quant-ph/0610161.
- D. McNulty, S. Weigert, "Mutually Unbiased Bases in Composite
  Dimensions -- A Review", Quantum 10, 2051 (2026), arXiv:2410.23997.
- P. H. J. Lampio, P. R. J. Ostergard, F. Szollosi, "Orderly generation
  of Butson Hadamard matrices", Math. Comp. 89 (2020), arXiv:1707.02287.
- G. Zauner, "Quantendesigns: Grundzuge einer nichtkommutativen
  Designtheorie", PhD thesis, Universitat Wien (1999).

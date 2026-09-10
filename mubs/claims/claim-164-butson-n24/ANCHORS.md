# External anchors

Every stage of the claim-164 pipeline reproduces a value published
independently of this project. A verifier can check each anchor without
executing any code from this repository.

## 1. Hadamard enumeration vs. the published Butson classification

Lampio, Ostergard and Szollosi classified Butson Hadamard matrices
BH(n, q) up to monomial equivalence ("Orderly generation of Butson
Hadamard matrices", Math. Comp. 89 (2020); arXiv:1707.02287; electronic
tables at the Aalto Butson wiki, `wiki.aalto.fi/spaces/Butson/`, page
"Matrices up to monomial equivalence"). Their table for n = 6 ends at
q = 17.

Closing the orbits of our 7,584 dephased row sets (and of the analogous
enumerations at smaller orders) under the full monomial-equivalence
residual group (`code/bh6_full_monomial_classes.py`) gives:

| N | dephased row sets | classes (ours) | classes (published) |
|---|---:|---:|---:|
| 3 | 12 | 1 | 1 |
| 4 | 72 | 1 | 1 |
| 6 | 312 | 4 | 4 |
| 8 | 432 | 3 | 3 |
| 12 | 2,184 | 11 | 11 |
| 24 | 7,584 | **25** | (none published) |

Agreement at every order with a published value. The count of 25 at
N = 24 (14 up to transpose) is new; the representatives are in `data/`.
The sub-Butson embedding is consistent: of the 25 classes, exactly 4
contain a BH(6, 6) representative, 3 a BH(6, 8), and 11 a BH(6, 12),
matching the published counts again, and 12 classes are new at root
order 24.

## 2. Completion census vs. Bengtsson et al. 2007

Bengtsson, Bruzda, Ericsson, Larsson, Tadej and Zyczkowski ("Mutually
unbiased bases and Hadamard matrices of order six", J. Math. Phys. 48,
052106 (2007); quant-ph/0610161), section 5, numerically enumerated all
order-6 Hadamard matrices with 24th-root entries, the vectors unbiased
to each, the third bases these admit, and whether any two third bases
form four MUBs. Their Fig. 2 census lists the equivalence classes that
admit at least one third basis, with candidate counts: F(0,0): 4,
F(1/6,0): 1, F^T(1/6,0): 1, F(1/6,1/12): 4, D(1/8): 4.

Our exact computation finds exactly five completion-bearing classes,
matching class for class:

| class profile (size, U rays, bases) | visible at sub-order | identification |
|---|---|---|
| 60, 12, 4 | 6, 12 | F(0,0) (the Fourier matrix) |
| 120, 6, 1 | 6, 12 | F(1/6, 0) |
| 120, 6, 1 | 6, 12 | F^T(1/6, 0) |
| 180, 24, 4 | 12 | F(1/6, 1/12) |
| 180, 24, 4 | 8 | D(1/8) (enphased Dita matrix) |

plus one class (72 members, 40 rays, 0 bases) with unbiased vectors but
no third basis, which their figure correctly omits. Their conclusion
(only triplets, no quadruple) is the numerical predecessor of this
claim; the present computation is exact and certified where theirs was
floating-point, and the two censuses agreeing class for class is a
strong mutual consistency check.

## 3. The published theorem record

McNulty and Weigert's review ("Mutually Unbiased Bases in Composite
Dimensions -- A Review", Quantum 10, 2051 (2026); arXiv:2410.23997)
states the exclusion as Theorem 7.15 only for BH(6, 12), citing the
same 2007 paper. The theorem record therefore stops at root order 12;
this claim extends it to 24 with a certified proof. Restricting our
pipeline to N = 12 reproduces the setting of Theorem 7.15 exactly
(2,184 dephased matrices, completion census {0/0: 1,884, 6/1: 240,
12/4: 60}, 480 completion bases, no unbiased pair).

## 4. Float margins (for re-implementers)

Anyone re-implementing the enumeration in floating point should know
the exact margins, computed in `code/butson_n24_referee_check.py` over
all 24^5 dephased rays: the smallest distance of a non-zero-sum ray
from zero-sum is 0.035, the smallest distance of a non-unbiased ray
from unbiasedness is 1.2e-3, and no ray falls in an ambiguous band at
double precision. Exact arithmetic (integer vectors modulo
Phi_24(x) = x^8 - x^4 + 1) avoids the question entirely and is what all
three certified implementations use.

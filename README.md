# Centaur Researcher

This repository holds notes for the "Centaur Researcher" project.

## Claim-164 (August 21, 2026)

Statement of theorem:

> There is no set of four mutually unbiased bases in C^6 of the form {I, H1/√6, H2/√6, H3/√6} in which every entry of H1, H2, H3 is a 24th root of unity. Equivalently: no MU quadruple in dimension six can be built from Butson Hadamard matrices in BH(6, 24).

Proof shape, in one sentence:
 
>  By exact arithmetic in Z[ζ24], the 7,584 dephased BH(6,24) matrices admit exactly 1,920 completion bases (candidate third bases unbiased to I and H1), and every one of the 2,520 candidate (H2, H3) pairs fails cross-unbiasedness, so the exhaustion closes.

Scope notes that must accompany any public statement:

> - Since ζm embeds in ζ24 for every m | 24, the theorem covers all root orders dividing 24; combined with your earlier exhaustion it covers every root order N ≤ 12 and N = 24. It says nothing about other root orders, mixed algebraic entries, or the continuous Hadamard families, and so does not resolve Zauner's N(6) = 3.
> - The stratum is not vacuous: it contains genuine MU triplets (including the exceptional 24th-root triplet), all 1,920 of which are shown not to extend.
> - Correct novelty framing: first exact, certified, independently cross-verified proof, confirming the uncertified numerical search of Bengtsson et al. 2007 (quant-ph/0610161, §5) and adding the first classification of BH(6,24) (25 monomial-equivalence classes, 14 up to transpose).


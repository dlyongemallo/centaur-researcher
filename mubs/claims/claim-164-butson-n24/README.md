# Claim 164: no MUB quadruple from Butson Hadamard matrices of root order 24

**Status: internally certified, awaiting external verification.**

## Statement

> **Theorem.** There is no set of four mutually unbiased bases of C^6 of the
> form {I, H1/sqrt(6), H2/sqrt(6), H3/sqrt(6)} in which every entry of H1,
> H2, and H3 is a 24th root of unity. Equivalently, no MU quadruple in
> dimension six can be built from Butson Hadamard matrices in BH(6, 24).

Since zeta_m embeds in zeta_24 for every m dividing 24, the theorem covers
all entries drawn from roots of unity of order dividing 24, in any mixture.
A corollary extends the exclusion to every single root order N <= 12 (for
N in {5, 7, 9, 10, 11} the relevant ray sets are empty, a seconds-long
computation included here; the rest divide 24).

**Scope.** The theorem says nothing about other root orders, mixed
algebraic entries outside mu_24, or the continuous families of complex
Hadamard matrices, and so does not resolve Zauner's conjecture that no
four MUBs exist in dimension six. The stratum is not vacuous: it contains
genuine MU triplets (including the exceptional 24th-root triplet), all
1,920 of which are shown not to extend.

**Novelty.** This is the first exact-arithmetic, certified, independently
cross-verified proof at root order 24. It confirms the uncertified
numerical search of Bengtsson et al. 2007 (quant-ph/0610161, section 5)
and adds the first classification of BH(6, 24): 25 monomial-equivalence
classes, 14 up to transpose (`data/`).

## Proof shape

By exact arithmetic in Z[zeta_24] (reduced modulo the cyclotomic
polynomial Phi_24(x) = x^8 - x^4 + 1), the 7,584 dephased BH(6, 24)
matrices are enumerated completely; they admit exactly 1,920 completion
bases (candidate third bases unbiased to I and to the matrix), and every
one of the 2,520 candidate (H2, H3) base pairs fails cross-unbiasedness.
The reduction of the theorem to this finite assertion is proved in
`note/`; the finite assertion is checked by the code in `code/`.

## Package contents

- `note/` -- the reduction note: the mathematics that equates the theorem
  with the finite computation. Read this first.
- `code/` -- the certified programs, byte-identical to the internally
  certified artefacts (see `SHA256SUMS`):
  - `butson_n24_independent.py` -- **primary reference implementation**
    (self-contained with numpy; ~2.5 min);
  - `butson_n24_search.py` -- the original search (~4.5 min);
  - `butson_n24_referee_check.py` -- adversarial referee harness with an
    independently re-derived reduction table and float-margin analysis
    (~6 min, three stages);
  - `butson_bh6_equivalence_classes.py`, `bh6_full_monomial_classes.py`,
    `make_class_representatives.py` -- equivalence-class census and the
    BH(6, 24) classification;
  - `butson_exhaust.py` -- the N <= 12 corollary (~25 s);
  - `mub_common.py`, `butson_n24_vanishing_sum.py` -- shared exact
    arithmetic.
- `logs/` -- the committed deterministic result logs (JSON).
- `data/` -- the 25 BH(6, 24) class representatives with completion
  census and sub-order visibility.
- `ANCHORS.md` -- the crosswalk to independently published data: every
  pipeline stage checks against a value in the literature.
- `SHA256SUMS` -- hashes of all code, log, and data files.

## Verify in ten minutes

Requires Python 3 with numpy. From this directory:

```sh
sha256sum -c SHA256SUMS

# Primary check: independent enumerator, full recomputation (~2.5 min).
# Regenerates logs/butson_n24_independent.json; it must differ from the
# committed version only in the "seconds" field.
python3 code/butson_n24_independent.py
git diff logs/butson_n24_independent.json

# Spot-check the committed search log and exact ray arithmetic (~15 s).
python3 code/butson_n24_search.py --verify
```

The expected census, at every stage: 8,350 zero-sum rays; 19,320
unbiased rays; 7,584 dephased Hadamards; completion census
{0 rays: 6,852, 6 rays/1 basis: 240, 12/4: 60, 24/4: 360, 40 rays/0
bases: 72}; 1,920 completion bases; 2,520 candidate pairs; 0 mutually
unbiased pairs.

Deeper verification, in increasing order of independence:

```sh
# Full original search (~4.5 min; writes to a historical path).
mkdir -p problems/mub-quad/logs
python3 code/butson_n24_search.py
diff <(python3 -c "import json;d=json.load(open('problems/mub-quad/logs/butson_n24_numpy_search.json'));d.pop('seconds');print(d)") \
     <(python3 -c "import json;d=json.load(open('logs/butson_n24_numpy_search.json'));d.pop('seconds');print(d)")

# Adversarial referee harness (~6 min; state under /tmp).
python3 code/butson_n24_referee_check.py --stage 1
python3 code/butson_n24_referee_check.py --stage 2
python3 code/butson_n24_referee_check.py --stage 3

# Classification census and the N <= 12 corollary.
python3 code/butson_bh6_equivalence_classes.py verify --log logs/butson_bh6_equivalence_classes.json
python3 code/butson_exhaust.py
```

The strongest verification needs none of this code: re-implement the
five-stage algorithm from the contract in `note/` (about a hundred lines
in any language with exact arithmetic) and compare your census against
the table above and the anchors in `ANCHORS.md`. Refutations are as
welcome as confirmations; see `../../../VERIFYING.md`.

## Provenance

Authored and cross-verified by three independent AI implementations
(Google, Anthropic, and OpenAI model families, 2026-08-21) under the
protocol described in `../../../PROVENANCE.md`; all three re-run cold
from a fresh clone by the human author on 2026-08-24 with byte-identical
results. The literature reconciliation (`ANCHORS.md`) was completed on
2026-08-24.

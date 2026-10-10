# Adversarial verification of the Codex reduction walk for CLAIM-164

- Date: 2026-10-06.
- Model: Claude Fable 5.1 (`claude-fable-5-1`).
- Review mode: adversarial verification of Codex's walk (`REDUCTION_WALK_codex.md`, 2026-10-05); read-only. No repository file other than this report was modified. All computations were run from a scratchpad against the committed code; `sha256sum -c SHA256SUMS` passes for every listed file.
- Files reviewed:
  - `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/note/REDUCTION_WALK_codex.md`
  - `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/note/butson_n24_note.tex`
  - `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/code/butson_n24_independent.py`
- Also consulted (supporting, not part of the reduction): `README.md`, `ANCHORS.md`, `SHA256SUMS`, `logs/*.json`, `data/bh6_24_class_representatives.txt`, `code/butson_n24_referee_check.py` (my 2026-08-21 harness), `code/butson_n24_search.py`, `code/butson_exhaust.py`, `code/mub_common.py`, `code/make_class_representatives.py`, `code/bh6_full_monomial_classes.py`, `code/butson_bh6_equivalence_classes.py`. Not read: the untracked `REDUCTION_WALK_agy.md`, to keep this pass independent of it.
- Disclosure: I authored the referee harness and signed off on CLAIM-164 at T2c. Finding B1 records what that sign-off did and did not cover.

## Summary

Codex's one concrete code finding (S2-1) is correct: `butson_n24_independent.py:114` tests completion rays against the rows of H_1, so the enumerator computes U(H_1^T) in the note's notation. I confirm this by line-by-line reading and by direct evaluation. I sharpen it in two directions Codex did not take: (i) the slip is provably a no-op for every aggregate the program reports, so Codex's "same digest after the fix" is a theorem rather than an observation; (ii) it is not a no-op per matrix: 720 of the 7,584 matrices (the F(1/6,1/12)/F^T and D(1/8)/D^T class pairs) swap profiles 24/4 and 0/0, and the external anchor (Bengtsson et al. Fig. 2) distinguishes the two conventions and sides with the note. The other two programs (search, referee harness) already use the note's convention, so the finite predicate as written is verified by two of the three programs; Codex's "the demanded object-by-object correspondence fails" is true of one file, not of the package. Codex's finding D (the abstract's "equivalently" overstates) is refuted: the standard basis is unbiased to every Butson basis, so Codex's scenario is already covered by the theorem; one clause would remove a residual misreading. The remaining findings are wording or hygiene. Verdict: nothing keeps the claim frozen; move toward ratification once line 114 is fixed and the artefacts regenerated.

Scorecard: S1-1 PARTIALLY CONFIRMED; S1-2 DOWNGRADED; S1-3 DOWNGRADED; S1-4 PARTIALLY CONFIRMED; S2-1 CONFIRMED (and sharpened); S2-2 CONFIRMED (minor); S2-3 CONFIRMED (minor); ordered/unordered DOWNGRADED; A, C, E agree (foreclosed); B agree modulo S1-1; D REFUTED; F agree with conclusion.

## Phase A: Codex's findings

Each entry: location, verdict, justification.

### S1-1. Row sorting omits the common action

- Note 149-161. PARTIALLY CONFIRMED (wording, not a gap).
- The proof's phrase "left multiplication by a permutation matrix, which fixes the standard basis and every MU relation" is the common-unitary argument: P applied to all four bases fixes I as a set and preserves every MU relation. What the enumeration needs is that the predicate "H_1 extends to a quadruple" is invariant under row permutations of H_1, which is exactly this. One clause ("applied to all four bases") fixes the lemma statement. Codex's extra step "re-dephase the columns of H_2, H_3 afterwards" is unnecessary: P fixes coordinate 0, so dephased columns stay dephased.

### S1-2. Cyclotomic-coordinate uniqueness asserted without proof

- Note 110-116. DOWNGRADED.
- "Z[omega] = Z[x]/Phi_24(x)" with Phi_24 monic of degree 8 is the statement that 1, x, ..., x^7 is a Z-basis; coefficientwise equality follows. This is textbook and appropriate for a note at this level; a citation is optional. I re-derived the table: the coefficients of x^k mod Phi_24 lie in {-1, 0, 1} for all k < 24, the code's recurrence x^8 = x^4 - 1 is correct, and the code's M-test (coefficient 0 equals 6, the other seven vanish) is the right test because |S|^2 is a rational integer.

### S1-3. Two linear-algebra facts uncited

- Note 119-122. DOWNGRADED.
- Both facts are stated with their one-line reasons already (|entry|^2 = 1/6; "both say H/sqrt6 is unitary"). Nothing to add.

### S1-4. Finite truth delegated to unnamed artefacts

- Note 218-256. PARTIALLY CONFIRMED.
- The package names and hashes the artefacts (`README.md`, `SHA256SUMS`, which passes), but the note's text does not. For a note meant to travel without the repository, name the three programs and pin a revision. Editorial.

### S2-1. The checker computes U(H_1^T)

- Code 110-116 against note 164-189. CONFIRMED.
- Line-by-line reading:
  - Line 110: `rows[i]` is row i of H_1 as an exponent vector, so `rows[i][j]` is entry (i, j).
  - Line 111: `cols = list(zip(*rows))` gives `cols[j] = (rows[0][j], ..., rows[5][j])`, i.e., column j; hence `cols[j][k] = rows[k][j]`.
  - Lines 113-116: for r in M and k = 1..5 the vector tested is `(r[j] - cols[j][k])_j = (r[j] - rows[k][j])_j`, i.e., r minus row k of H_1.
  - Note 165-169 defines U(H_1) by c - h^(i) in M for the columns h^(i) = `cols[i]` = `(cols[i][j])_j`; the matching code is `cols[k][j]`.
  - Direct evaluation (scratchpad, importing the enumerator's own functions): on the first 50 matrices the literal expression of lines 113-116 equals the set {c in M : c - r in M for all r in R} exactly. The enumerator therefore computes U(H_R^T). Codex's one-token fix is right.
- Why every reported aggregate is provably unaffected. Codex asserted this from one recomputation; here is the argument. Let f(H) = (|U(H)|, number of 6-cliques, pair verdict). f is invariant under column permutations of H (U depends only on the column set) and under row permutations P fixing row 0 (U(PH) = P.U(H), and Z, M are invariant under coordinate permutations fixing coordinate 0, so cliques and M-tests transport). Let 𝓗 be the set of all 7,584 x 120 dephased BH(6,24) matrices. Transposition T is a bijection of 𝓗, and f is constant on row-orbits (the sets R) and on column-orbits. Hence sum_R f(H_R^T) = (1/120) sum_{H in 𝓗} f(TH) = (1/120) sum_{H in 𝓗} f(H) = sum_R f(H_R). Census, total completions, `max_u`, `max_completions`, rejected-pair count, verdict and the stable digest are identical under the two conventions. Codex's observation that the digest `2eac89c3...` survives the fix is a theorem, not evidence.
- Why it is nevertheless not a no-op (new; Codex did not check). Computed today over all 7,584 matrices, the joint census under the (note, code) conventions is: 0/0 | 0/0: 6,492; 24/4 | 0/0: 360; 0/0 | 24/4: 360; 12/4 | 12/4: 60; 6/1 | 6/1: 240; 40/0 | 40/0: 72. So 720 matrices change profile and 1,090 have U(H_R) != U(H_R^T) as sets. The 720 are classes 6, 17, 20, 21 of `data/bh6_24_class_representatives.txt`: the F(1/6,1/12) / F^T(1/6,1/12) pair (21 / 6) and the D(1/8) / D^T(1/8) pair (20 / 17), with 24/4 on one member and 0/0 on the other. The data file already records this (class 20: U=24, bases=4, transpose_class=17; class 17: U=0), so the per-matrix census is manifestly not transpose-invariant.
- Why no committed artefact is wrong. The independent enumerator reports aggregates only. The per-class data file is generated by `make_class_representatives.py` through `be.completion_rays`, which tests columns (`columns = list(zip(*h_rows))`, note convention). My referee harness tests columns at lines 239-246 (`cols[1:, None, :]`), cross-checks against `be.completion_rays` (lines 250-254), and verifies each of the 1,920 bases exactly against its H_1 with `mc.exact_mu_pair(h1, b, D)` (line 319), i.e., |(H_1^* H_2)_kl|^2 = 6, the column-basis relation. The search uses `be.completion_rays`. So the note's predicate is verified exactly as written by two programs and in transposed form by the third. Codex's sentence "the demanded object-by-object correspondence fails" is true of `butson_n24_independent.py` and false of the package.
- Which convention is externally right (new). Constructing F_6^(2)(a, b) from the Tadej-Zyczkowski pattern (phases (0, a, b, 0, a, b) on rows 1, 3, 5; a = 1/6 gives exponent 4, b = 1/12 gives exponent 2 in Z/24) and closing orbits with the package's own monomial moves (`bh6_full_monomial_classes.neighbours`) gives: F(0,0) in class 0 (orbit 60); F(1/6,0) in class 12 and F^T(1/6,0) in class 4 (both 6/1); F(1/6,1/12) in class 21 (24/4) and F^T(1/6,1/12) in class 6 (0/0). Bengtsson et al. Fig. 2 lists F(1/6,1/12) and not its transpose, so the column convention is the one anchored to the literature; under the row convention the attribution reverses and contradicts Fig. 2. The fix therefore belongs in the code, not in the note. Caveat: this relies on my recollection of which of the two transposed two-parameter Fourier families the catalogue calls F_6^(2); check against Tadej-Zyczkowski 2006 before quoting.
- Codex's rescue argument is correct in content for the existence verdict. The equality of the full census additionally needs the observation that f is a class function on both row- and column-orbits, given above. After the fix, `SHA256SUMS` and the README sentence "byte-identical to the internally certified artefacts" must be updated; the log's deterministic digest does not change.

### S2-2. `--verify` recomputes nothing

- Code 152-166. CONFIRMED (minor).
- Accurate description. The README never directs a verifier to this flag (it directs to the full run plus `git diff` of the log), and the search's analogous flag is labelled a spot-check, so no reader of the package is misled. Renaming to `--verify-log` is sensible.

### S2-3. Default run does not assert the later counts

- Code 99-103, 134-150. CONFIRMED (minor).
- `main()` asserts |Z|, |M| and the Hadamard count only. The referee harness asserts every stage count against `COMMITTED` (lines 70-78, 259, 342-345), so the package has asserting code; adding three asserts here costs nothing.

### Ordered versus unordered pairs

- Note 237-243; code 123-124. DOWNGRADED.
- The note writes 420 x C(4,2) = 2,520, which already says unordered. Adding the word is harmless.

### A (missed stratum), C (pair counting), E (normalisation constant)

- Agree: foreclosed, for the reasons Codex gives. For E I checked that the code's M-test is |S|^2 = 6 in integer coefficients and its Z-test is S = 0, matching the note's |<u/sqrt6, w/sqrt6>|^2 = |S(w - u)|^2 / 36.

### B (hidden dephasing freedom)

- Agree: foreclosed modulo the S1-1 wording.

### D. The "equivalently" clause broadens the theorem

- Note 28-32. REFUTED as argued; one residual editorial point.
- Codex's scenario is four pairwise unbiased bases K_0, ..., K_3 with K_i in BH(6,24) and I not among them, and the worry is that K_0^* K_i / sqrt6 need not have root-of-unity entries. But I is automatically unbiased to every K_i / sqrt6, since every entry has modulus 1/sqrt6 (the note says so at 119-120), so {I, K_1/sqrt6, K_2/sqrt6, K_3/sqrt6} is already a quadruple of the theorem's form. No change of frame is needed, and Codex's "Status: Open" is wrong. The theorem implies the broad statement and in fact a stronger one: no three pairwise unbiased bases given by BH(6,24) matrices exist in any common frame, since any such triple extends by I. The one reading not covered is "each H_i Hadamard-equivalent to a BH(6,24) matrix", because Hadamard equivalence allows a left diagonal with arbitrary phases; a clause such as "with all entries of the three transition matrices in mu_24" precludes it. In my referee role: the abstract does not overstate the theorem. Codex's proposed fix ("define 'built from' to mean the identity-containing form") is more restrictive than the theorem warrants; recommend instead adding the free-identity corollary as a remark and the clarifying clause.

### F. The transpose slip could hide a matrix-specific completion

- Agree with the conclusion (universal verdict foreclosed); see S2-1 for the proof. Per-matrix data do move, so the concern was legitimate to raise.

## Phase B: independent findings

- B1. Harness coverage and what my August sign-off missed. The 2026-08-21 harness verifies the note's predicate in the column convention for all 7,584 matrices, exactly checks all 1,920 bases (orthonormality and unbiasedness to H_1) and all 2,520 pairs in a second exact representation, re-derives the reduction table, and cross-checks U(H_1) against the search's `completion_rays`. It never compared against the independent enumerator's U sets (the enumerator exposes none), so the convention mismatch was invisible to it, and I did not catch it at T2c. The harness obviates Codex's strongest phrasing of S2-1 but not the one-token defect itself.

- B2. Per-matrix non-invariance (numbers in S2-1). Consequence for the choice of fix: correct the index rather than document the transpose, because the alternative leaves the "primary reference implementation" computing per-matrix data that disagree with `data/` on 720 matrices.

- B3. Anchor identification for the transpose-asymmetric classes. `ANCHORS.md` section 2 and note section 5 identify class 20 as D(1/8) and class 21 as F(1/6,1/12). Classes 20/17 share (size 180, sub-order 8) and 21/6 share (size 180, sub-order 12); the package contains no code matching catalogue matrices, so the choice of member within each pair could only have come from the completion count, which is what the anchor is meant to test. Today's orbit computation resolves the F pair independently (F(1/6,1/12) in class 21, with the catalogue-pattern caveat). The D pair is consistent with Fig. 2 but not independently constructed; a one-time explicit construction of D_6^(1)(1/8) would close it. Supporting material, not part of the reduction.

- B4. "1,920 MU triplets" (note 85-87; README). 1,920 counts completion bases, i.e., normal-form pairs (H_1, H_2) with H_1 row-sorted and dephased and H_2 a 6-clique. An unordered triplet {I, A, B} is reached once from A and once from B, so the number of unordered triplets is about half. Say "1,920 completion bases (MU triplets with a distinguished first Hadamard)".

- B5. Lemma 3 statement, note 171-173: "representable with columns in BH(6,24)" should read "representable by a matrix in BH(6,24)".

- B6. Corollary 2 sanity check (floating point, not a certificate): min |S| over dephased 6-tuples is 0.38, 0.25, 0.088 at N = 5, 7, 11 (no vanishing sums, as Lam-Leung predicts), and min ||S|^2 - 6| is 0.087, 0.090 at N = 9, 10 (the analogue of M is empty). Consistent with the note.

- B7. Code review of `butson_n24_independent.py` beyond line 114: reduction table, ray classification (numpy promotes int16 sums to int64, and the sums are bounded in any case), packed keys (24^6 < 2^28), upper-triangular adjacency with increasing-order clique search, empty-U handling, and the pair loop are all correct. Nothing else found.

- B8. Process consequence of the fix: changing line 114 changes the code hash (`SHA256SUMS`, README provenance sentence) but not the log's deterministic digest; the note's "three independently written programs ... all agree" then holds in the note's convention for all three.

## Phase C: verdict

Nothing in Codex's walk or in this pass touches the theorem, the reduction, or the finite assertion. The predicate of Proposition 1 is verified exactly as stated by two programs, one of them with exact per-base and per-pair checks in a second representation, and in a provably equivalent transposed form by the third.

Required before ratification: fix `butson_n24_independent.py:114` to `cols[k][j]`; regenerate `logs/butson_n24_independent.json` (digest unchanged) and `SHA256SUMS`; update the README provenance sentence; add the S1-1 clause to Lemma 2.

Recommended: name the artefacts in the note (S1-4); add the free-identity remark and a clarifying clause to the abstract (D); fix the "1,920 triplets" wording (B4) and the Lemma 3 wording (B5); rename `--verify` to `--verify-log` and add the three asserts (S2-2, S2-3); construct D(1/8) explicitly for the anchor (B3).

CLAIM-164 should move toward ratification once the required items are done. This pass finds nothing that keeps it frozen beyond the "awaiting ratification" state.

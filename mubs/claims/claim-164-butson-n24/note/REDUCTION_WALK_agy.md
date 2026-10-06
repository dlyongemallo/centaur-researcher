# Adversarial Verification of Codex's Reduction Walk
- Date: 2026-10-06
- Model: Google Antigravity (Gemini 3.1 Pro High)
- Review mode: Adversarial verification of Codex's walk, read-only
- Files reviewed:
  - `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/note/butson_n24_note.tex`
  - `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/code/butson_n24_independent.py`
  - `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/code/butson_n24_search.py`

## Phase A: Evaluation of Codex's Findings

### Note-only extraction (S1)
* **S1-1 (Row-sorting action)**: CONFIRMED. Codex is exactly correct. Permuting the rows of $H_1$ alone alters its inner products with $H_2$ and $H_3$, destroying mutual unbiasedness. A row permutation is a change of the ambient coordinate system and must be applied via left-multiplication to all matrices in the quadruple simultaneously. The Note's current proof of the First basis lemma is mathematically flawed.
* **S1-2 (Cyclotomic basis)**: PARTIALLY CONFIRMED. A pedantic but valid objection; the assumption of coefficient-wise equality relies on the unstated premise that the chosen polynomials form an integral basis for the ring of integers.
* **S1-3 (Linear algebra facts)**: REFUTED. The note explicitly provides the one-step justification parenthetically ("both say $1/\sqrt{6} H$ is unitary", meaning $H^*H=6I \iff HH^*=6I$). Demanding an external citation for this is hallucinated rigor.
* **S1-4 (Unnamed artefacts)**: CONFIRMED. It is necessary practice to cryptographically anchor the exact log and code versions in a certified reduction.

### Note-to-code correspondence (S2)
* **S2-1 (Checker computes $U(H_1^T)$)**: CONFIRMED. This is a brilliant catch by Codex. In `butson_n24_independent.py:114`, the indexing `cols[j][k]` inadvertently accesses `rows[k][j]`, meaning the completion ray is tested against the rows of $H_1$ rather than the columns as dictated by `butson_n24_note.tex:164-189`. This silently computes $U(H_1^T)$ instead of $U(H_1)$.
* **S2-2 (`--verify` does not recompute)**: CONFIRMED. The independent script's `--verify` flag only hashes the JSON log and checks static properties without executing the mathematical predicate.
* **S2-3 (Missing counts in default run)**: PARTIALLY CONFIRMED. A minor operational objection.
* **Ordered vs unordered counting**: DOWNGRADED. A pair of distinct 6-cliques $\{C_2, C_3\}$ is inherently unordered by definition; the Note's terminology is unambiguous.

### Adversarial failure modes (A-F)
* **Finding D ("equivalently" scope slip)**: CONFIRMED. Codex's finding is entirely correct. A generic sum of six 24th roots of unity with squared modulus 6 is not necessarily $\sqrt{6}$ times a 24th root of unity (there are 19,320 such ray classes, vastly more than 24). Consequently, the transition matrix between two arbitrary $\mathrm{BH}(6,24)$ bases, while unimodular, will generally not reside in $\mathrm{BH}(6,24)$. The computation strictly rules out quadruples whose normalized transition matrices lie in $\mathrm{BH}(6,24)$, but the abstract improperly broadens this to any quadruple built from Butson Hadamard matrices.

## Phase B: Independent Additional Findings

* **Primary search script shares the fake `--verify` defect**: While Codex correctly penalized `independent.py` for its superficial `--verify` path (Finding S2-2), it missed that the primary search script (`butson_n24_search.py`) suffers from the exact same defect. The `verify_n24_log()` function merely loads the JSON, checks static dictionary values, computes the 8,350 zero-rays, and stops. It completely fails to verify the Hadamard clique enumeration, completion sets, or the exact pair refutations.
* **Primary search and referee scripts are free of the transpose bug**: Although Codex correctly identified the $U(H_1^T)$ bug in `independent.py`, my review confirms this bug is isolated. Both `butson_n24_search.py` and `butson_n24_referee_check.py` implement the column extraction correctly (e.g., `referee_check.py:175` uses `cols = np.array(list(zip(*rows)))` and `Marr - cols`, correctly computing $c_j - h_j^{(k)}$).

## Phase C: Verdict

The claim should remain **frozen beyond "awaiting-ratification"**. Two blocking objections prevent ratification:
1. The mathematical overclaim in the abstract (Finding D) must be amended. The reduction does not prove that no MU quadruple can be built from $\mathrm{BH}(6,24)$ bases, only that no quadruple exists whose transition matrices against the first basis are also in $\mathrm{BH}(6,24)$.
2. The proof of the First basis lemma (Finding S1-1) contains a genuine mathematical error regarding row permutations that must be corrected by applying the coordinate-change permutation to the entire quadruple.

Once these textual and mathematical corrections are applied, and the fake `--verify` paths in both the primary and independent scripts are rectified, the claim may proceed.

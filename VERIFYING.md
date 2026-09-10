# Verifying claims in this repository

## What a claim package contains

Each claim directory (e.g., `mubs/claims/claim-164-butson-n24/`) is
self-contained:

- `README.md` -- the exact statement, its scope, and a quickstart.
- `note/` -- a short mathematical note reducing the statement to a
  finite, machine-checkable assertion. This is the part to *read*; no
  amount of re-running code substitutes for checking it.
- `code/` and `logs/` -- programs that check the finite assertion, and
  the deterministic logs they are expected to reproduce. `SHA256SUMS`
  ties the published files to the internally certified artefacts.
- `ANCHORS.md` -- points of contact with independently published data,
  so each pipeline stage can be checked against the literature without
  trusting our code.

## Levels of verification

In increasing order of strength:

1. **Re-run**: execute the packaged code and confirm it reproduces the
   committed logs (differences confined to wall-clock fields). Minutes.
2. **Cross-check**: verify the anchors against the cited literature,
   and run the independent implementations against each other.
3. **Re-implement**: write your own program from the algorithm contract
   in the note, without reading our code, and compare censuses. This is
   the verification we most value, and the notes are written to make it
   an afternoon's work.
4. **Audit the mathematics**: check the reduction lemmas in the note.
   For nonexistence claims, this step is essential: the code only
   checks the finite assertion, and the note is what proves the finite
   assertion equals the theorem.

## Reporting

Please report outcomes as GitHub issues on this repository:

- **Confirmation**: state which level(s) you performed, your
  environment, and the census/hashes you obtained. Confirmations are
  credited in the claim's README unless you prefer otherwise.
- **Refutation or gap**: a counterexample, a census mismatch, or an
  error in a reduction lemma. Refutations are as valuable to us as
  confirmations and will be published just as prominently, with credit.
  Where possible include the object (e.g., an explicit matrix) or the
  step of the note you dispute.

Questions and partial results are welcome in issues as well.

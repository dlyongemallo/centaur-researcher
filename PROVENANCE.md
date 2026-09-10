# Provenance

The results in this repository are produced by a "centaur" workflow: an
autonomous research loop of large language models from three vendors
(Anthropic, Google, and OpenAI model families) working over the human
author's private research repositories, directed and audited by the
human author (David Yonge-Mallo). This page discloses how a claim earns
publication here; the disclosure is part of the evidence.

## The internal protocol

- **Tiered claims ledger.** Every result is recorded with committed,
  re-runnable evidence. No session may build on a claim beyond its
  verification tier.
- **Cross-family refereeing.** A claim is promoted only after a session
  of a *different* model family re-runs its certificate cold, from a
  fresh clone, and reproduces the committed logs byte for byte (modulo
  wall-clock fields). Referees are prompted adversarially and log
  objections; rejections are sticky.
- **Milestone freezing.** Result-level claims are additionally frozen
  until the human author personally re-runs the certificates and
  ratifies. Nothing internal may build on an unratified milestone.
- **Exactness discipline.** Floating-point output never counts as a
  theorem; certificates are exact (integer/rational arithmetic,
  cyclotomic reduction, or outward-rounded interval arithmetic).

## Claim 164 specifically

- 2026-08-21: search authored by a Google-family model (exact
  Z[zeta_24] arithmetic with vectorized clique assembly). The same day,
  an Anthropic-family referee session wrote an independent verification
  harness (hand-derived reduction table, independent enumeration and
  margin analysis over all 24^5 rays, exact verification of all bases
  and refuted pairs in a second exact representation) and reproduced
  the search log cold; an OpenAI-family session then wrote a third,
  code-independent enumerator reproducing the full census.
- 2026-08-24: the human author re-ran all three implementations from a
  fresh clone; all logs reproduced byte-identically modulo wall-clock
  fields. The literature reconciliation in the claim's `ANCHORS.md`
  (including the corrective finding that Bengtsson et al. 2007 already
  searched 24th roots numerically) was completed the same day.
- The code published in the claim's `code/` directory is byte-identical
  to the internally certified artefacts; `SHA256SUMS` in the claim
  directory pins the correspondence.

## What is deliberately not published

The internal ledger, task queues, dead-end maps, and unrelated research
threads are not part of any claim's evidence and are omitted to keep
packages reviewable. The full internal ledger for a published claim is
available to serious reviewers on request.

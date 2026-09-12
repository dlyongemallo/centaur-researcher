# Centaur Researcher

Selected results from a human-directed, AI-executed ("centaur")
research programme in quantum information and adjacent mathematics,
published here for independent verification or refutation by human
experts. Each claim is a self-contained evidence package: a precise
statement, a short proof note reducing it to a finite computation,
certified code and deterministic logs, and a crosswalk to independently
published data.

## Motivation

> Everyone's one-shotting major open problems with AI and I'm just sitting
> here working on minor toy problems by hand.

If you've followed AI developments, you'll know that the frontier labs have
been solving major open problems using LLMs. Like many researchers, I have a
lot of research notes on and attempts at solving various open problems.
As I am not in academia (or otherwise getting paid to work on this research),
there was very little chance that any of it would amount to anything...
until the age of LLMs dawned upon us.
I figured, why not set a bunch of agents loose on my research notes
and see what comes out of it? As I wanted to be a "centaur" rather than a
"reverse centaur"[^centaur], I set up a system where the LLMs would do the
work but humans-in-the-loop would still do the thinking and guiding.
You're looking at that project.

## Verification and Provenance

How verification works, and how to report a confirmation or a
refutation: see [VERIFYING.md](VERIFYING.md). How these results were
produced and cross-checked before publication: see
[PROVENANCE.md](PROVENANCE.md).

## Claims

| problem | claim | statement (short) | status |
|---|---|---|---|
| [Mutually unbiased bases in dimension six](mubs/) | [claim-164](mubs/claims/claim-164-butson-n24/) | No four mutually unbiased bases of C^6 can be built from Butson Hadamard matrices with 24th-root-of-unity entries. First certified proof at root order 24; includes the first classification of BH(6, 24) (25 equivalence classes). | internally certified, awaiting external verification |

Statuses: **internally certified** (independent implementations by
multiple AI model families reproduce the result exactly, ratified by
the human author) -> **externally confirmed** (at least one independent
human-directed reproduction, credited in the claim) or **refuted**
(a verified counterexample or error; refutations are welcome and will
be published just as prominently).

## Licensing

Code is released under the MIT License (see [LICENSE](LICENSE)); text,
notes, and data files under CC BY 4.0. Cite via
[CITATION.cff](CITATION.cff).

[^centaur]: Cory Doctorow, "The Reverse Centaur's Guide to Life After AI", Verso Books, 2026.

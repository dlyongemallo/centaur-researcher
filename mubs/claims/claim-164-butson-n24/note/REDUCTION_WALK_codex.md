# Adversarial Reduction Walk for Claim 164

> *Terminology note: this file is named `REDUCTION_WALK_codex.md` after the auditor model (OpenAI's Codex / gpt-5.6-sol), following the sibling convention `REDUCTION_WALK_<family>.md`. It is unrelated to the September 2026 `openai/math` preprint dump.*

- Date: 2026-10-05.
- Model: OpenAI Codex (GPT-5).
- Review mode: Independent, adversarial, and read-only.
- Note reviewed: `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/note/butson_n24_note.tex`.
- Checker reviewed: `/home/davinci/workspace/centaur-researcher-wt/main/mubs/claims/claim-164-butson-n24/code/butson_n24_independent.py`.

## Principal objections

1. `butson_n24_independent.py:110-116` tests completion rays against the rows of \(H_1\), not its columns as required by `butson_n24_note.tex:164-175`. It therefore computes \(U(H_1^T)\), not the stated \(U(H_1)\). Transposition symmetry rescues the global universal nonexistence result, but that extra argument is absent and the demanded object-by-object correspondence fails.
2. The “equivalently” language at `butson_n24_note.tex:28-32` is too broad unless “built from Butson Hadamard matrices” is explicitly defined to mean the displayed identity-containing form. Sending an arbitrary first Butson basis to \(I\) need not preserve the \(24\)th-root property of the other transition matrices.
3. The row-sorting justification at `butson_n24_note.tex:149-161` omits that the row permutation must be applied to every basis matrix, not only to \(H_1\).
4. The checker’s `--verify` path at `butson_n24_independent.py:152-166` verifies a stored JSON log rather than recomputing the finite predicate.

## Step 1: Note-only extraction

### The finite predicate

Let
\[
A=\{0\}\times(\mathbb Z/24)^5,
\qquad
S(a)=\sum_{j=0}^{5}\zeta_{24}^{a_j}.
\]
Define
\[
Z=\{a\in A:S(a)=0\},
\qquad
M=\{a\in A:S(a)\overline{S(a)}=6\}.
\]

Let \(\mathcal R\) be the explicitly finite collection
\[
\mathcal R=
\left\{
R\subseteq Z:
|R|=5,\ 
r'-r\in Z\text{ for all distinct }r,r'\in R
\right\}.
\]

For \(R\in\mathcal R\), form \(H_R\) with exponent rows consisting of the zero vector and the five elements of \(R\), in a fixed canonical order. Write its exponent columns as
\[
h^{(0)}=0,h^{(1)},\ldots,h^{(5)}.
\]
Define
\[
U_R=
\left\{
c\in M:
c-h^{(i)}\in M\text{ for }i=1,\ldots,5
\right\}.
\]
Define the completion bases
\[
\mathcal C_R=
\left\{
C\subseteq U_R:
|C|=6,\ 
c'-c\in Z\text{ for all distinct }c,c'\in C
\right\}.
\]

The exact finite assertion is
\[
\forall R\in\mathcal R,\ 
\forall C_2\ne C_3\in\mathcal C_R,\ 
\exists c\in C_2,\ c'\in C_3
\text{ such that }c'-c\notin M.
\]

Everything ranges over subsets of \(A\), whose size is exactly \(24^5\), so this is a self-contained predicate on an explicitly enumerable finite set.

Locations: `butson_n24_note.tex:91-117`, `butson_n24_note.tex:144-209`.

### Lemmas used by the reduction

The proof of the finite-reduction proposition at `butson_n24_note.tex:204-214` uses exactly these four lemmas:

1. Normal form, `butson_n24_note.tex:124-142`.
2. First basis, `butson_n24_note.tex:144-162`.
3. Second and third bases, `butson_n24_note.tex:164-189`.
4. Pair test, `butson_n24_note.tex:191-202`.

### S1-1: The row-sorting argument omits a common action

- Location: `butson_n24_note.tex:149-161`.
- What is missing: Permuting the rows of \(H_1\) alone does not generally preserve its mutual unbiasedness with fixed \(H_2\) and \(H_3\). The permutation must act on the left of every basis matrix simultaneously.
- Why it matters: The enumeration replaces the five noninitial rows by an unordered subset. The stated proof must justify that this quotient does not lose quadruples.
- Smallest fix: Replace the last sentence of the proof with an explicit statement that a permutation fixing the first coordinate is applied on the left to \(I,H_1,H_2,H_3\), after which columns of \(H_2,H_3\) may be re-dephased by right diagonal matrices.

### S1-2: Cyclotomic-coordinate uniqueness is asserted without proof or citation

- Location: `butson_n24_note.tex:110-116`.
- What is missing: The note states \(\Phi_{24}(x)=x^8-x^4+1\) and that ring elements are represented by integer vectors of length eight, but does not explain that \(1,\zeta_{24},\ldots,\zeta_{24}^7\) is an integral basis for the quotient representation and that equality is therefore coefficientwise.
- Why it matters: The claim that the program’s integer comparisons are exact rests on this representation.
- Smallest fix: Add a one-sentence derivation or a standard cyclotomic reference, and state that the monic degree-eight minimal polynomial gives unique representatives of degree below eight.

### S1-3: Two elementary linear-algebra facts are uncited

- Location: `butson_n24_note.tex:119-122`.
- What is missing: The note uses, without proof or citation, that a root-valued normalized column is unbiased to the standard basis and that a square matrix with orthogonal columns has orthogonal rows.
- Why it matters: These facts connect the row enumeration and implicit identity basis to the basis formulation.
- Smallest fix: Add one sentence: the first claim follows from every normalized entry having modulus \(1/\sqrt6\), and the second from \(H^*H=6I\) implying \(HH^*=6I\) for a square matrix.

### S1-4: The finite truth is delegated to unnamed artefacts

- Location: `butson_n24_note.tex:218-256`.
- What is missing: The numerical counts and all \(2{,}520\) pair refutations are asserted as computational results. The note mentions accompanying programs and manifests but gives no filenames, versions, commit identifier, or hashes.
- Why it matters: The finite assertion is the only non-elementary premise remaining after the reduction.
- Smallest fix: Identify the checker files and committed logs explicitly, and include either their SHA-256 values or a stable repository revision.

## Step 2: Note-to-code correspondence

### Crosswalk

| Required correspondence | Note | Code | Result |
|---|---|---|---|
| (a) Dephasing gauge | `butson_n24_note.tex:124-175` | `butson_n24_independent.py:38-55`, `101-118` | The finite rays have first exponent zero; \(H_1\) has a zero first row and first column; completion columns are dephased. The row-sorting proof needs the common-permutation clarification from S1-1. |
| (b) Basis of \(\mathbb Z[\zeta_{24}]\) | `butson_n24_note.tex:110-116` | `butson_n24_independent.py:19-32`, `47-52` | The code uses the power basis \(1,x,\ldots,x^7\) and the recurrence \(x^8=x^4-1\). Zero and norm-six tests are coefficientwise integer comparisons. |
| (c) Completion rays \(U(H_1)\) | `butson_n24_note.tex:164-189` | `butson_n24_independent.py:110-116` | Mismatch: the code uses rows where the note requires columns. |
| (d) Cross-unbiasedness | `butson_n24_note.tex:191-209` | `butson_n24_independent.py:122-133` | The code tests every cross pair until it finds a witness. Its difference direction is reversed, but \(M\) is invariant under negation, so this does not change the predicate. |
| (e) Meaning of quadruple | `butson_n24_note.tex:73-75`, `204-213` | `butson_n24_independent.py:109-133` | \(I\) is implicit, \(H_1\) is distinguished, and the pair \(\{H_2,H_3\}\) is unordered. This is sufficient for nonexistence of the set denoted by braces, but it is not an enumeration of labelled ordered quadruples. |

### S2-1: The checker computes \(U(H_1^T)\)

- Locations: `butson_n24_note.tex:165-175`; `butson_n24_independent.py:110-116`.
- What is wrong: After
  `cols = list(zip(*rows))`,
  the exponent vector of column \(k\) is
  \[
  (\texttt{cols[k][0]},\ldots,\texttt{cols[k][5]}).
  \]
  The checker instead forms
  \[
  (\texttt{cols[0][k]},\ldots,\texttt{cols[5][k]}),
  \]
  which is row \(k\) of the original matrix. Thus line 114 should use `cols[k][j]`, not `cols[j][k]`.
- Why it matters: The program does not have the required line-by-line correspondence with the formal \(U(H_1)\). Per-matrix completion data are assigned to the transpose rather than to the enumerated matrix.
- Why the global verdict survives: Transposition is a bijection of dephased \(\mathrm{BH}(6,24)\). Sorting the noninitial rows after transposition is a coordinate permutation fixing the dephased coordinate, and completion-pair existence is invariant under that common coordinate permutation. Therefore
  \[
  \forall H\,P(H^T)
  \quad\Longleftrightarrow\quad
  \forall H\,P(H).
  \]
  This rescue argument is not stated in the note or code.
- Further check performed: The unmodified checker was recomputed read-only and produced the reported census. An in-memory run replacing only `cols[j][k]` with `cols[k][j]` produced the same census, verdict, and stable digest:
  `2eac89c3d837af7162f3a5ae88f374362824f0fbd46382ff31b04bb8f8679506`.
- Smallest fix: Change `butson_n24_independent.py:114` to `cols[k][j]`, regenerate the log and manifests, and retain the corrected implementation. Alternatively, document and prove the transpose reindexing, but that leaves the requested object-by-object correspondence false.

### S2-2: `--verify` does not recompute the finite assertion

- Location: `butson_n24_independent.py:152-166`.
- What is wrong: `verify()` reads the previously generated JSON, checks expected fields, and recomputes a digest stored in that same JSON. It does not execute ray, Hadamard, completion, or pair enumeration.
- Why it matters: A user invoking `--verify` has verified a log’s internal consistency, not the mathematical predicate.
- Smallest fix: Rename the option to `--verify-log`, or make `--verify` recompute and compare the result with the committed log.

### S2-3: The default run does not assert all claimed checkpoint counts

- Locations: `butson_n24_note.tex:231-243`; `butson_n24_independent.py:99-103`, `134-150`.
- What is missing: The default computation asserts the ray and Hadamard counts, but it only prints the completion census, total completions, and rejected-pair count. The expected values for those later checkpoints occur only in the non-recomputing `verify()` path.
- Why it matters: The note calls the stage counts part of the claim and checkpoints for independent reimplementation.
- Smallest fix: Before writing the result, assert the expected completion census, total \(1{,}920\), and rejected-pair count \(2{,}520\).

### Ordered versus unordered counting

- Locations: `butson_n24_note.tex:204-209`, `237-243`; `butson_n24_independent.py:123-124`.
- Observation: The code takes `ia < ib`, so the \(2{,}520\) candidate pairs are unordered. An ordered \((H_2,H_3)\) count would be \(5{,}040\).
- Why it matters: This is not an undercount for nonexistence because cross-unbiasedness is symmetric. It should nevertheless be stated so that “candidate pair” is not mistaken for a count of labelled ordered quadruples.
- Smallest fix: Change the note to “\(2{,}520\) unordered candidate pairs.”

## Step 3: Adversarial failure modes

### A. A non-dephased or exceptional algebraic stratum is missed

- Locations: `butson_n24_note.tex:91-117`, `124-162`; `butson_n24_independent.py:38-55`, `101-103`.
- Status: Foreclosed for the narrow theorem.
- Reason: \(\mathrm{BH}(6,24)\) is finite. The Normal form lemma dephases \(H_1\) using only \(24\)th-root phases, and the First basis lemma characterizes every resulting matrix as a five-clique in \(Z\). The code enumerates all \(24^5\) dephased rays before enumerating all five-cliques. There is no continuous or measure-zero sub-stratum inside this finite set.
- Smallest further check: None beyond correcting or documenting the common row-permutation argument.

### B. Hidden dephasing freedom creates unenumerated quadruples

- Locations: `butson_n24_note.tex:124-142`, `164-189`; `butson_n24_independent.py:38-55`, `110-118`.
- Status: Foreclosed except for the row-sorting wording defect.
- Reason: Right diagonal phases put every column of \(H_2,H_3\) into exactly one dephased ray class, and column order is removed by six-element subsets. Common left diagonal phases preserve the identity basis as a set and keep entries in \(\mu_{24}\).
- Remaining issue: Sorting \(H_1\)’s rows must be described as a common left permutation of all four bases.
- Smallest fix: Amend the First basis proof as in S1-1.

### C. Candidate-pair overcounting or undercounting

- Locations: `butson_n24_note.tex:191-209`, `237-243`; `butson_n24_independent.py:117-133`.
- Status: Foreclosed for the existence predicate.
- Reason: Each completion basis is an unordered six-element set of dephased rays. The code considers each distinct unordered pair once. Swapping \(H_2\) and \(H_3\) cannot change whether all 36 cross-inner-products have the required modulus.
- Documentation objection: The note should call the \(2{,}520\) pairs unordered.

### D. The identity-basis assumption broadens the theorem improperly

- Locations: `butson_n24_note.tex:28-32`, `55-58`, `73-76`.
- Status: Open for the abstract’s broad “equivalently” formulation; closed only for the theorem exactly as displayed.
- Reason: Any MU quadruple can be sent to one containing \(I\) by a common unitary. However, if the original four basis matrices are \(K_0,K_1,K_2,K_3\in\mathrm{BH}(6,24)\), the resulting unnormalised transition matrix is
  \[
  H_i'=\frac{1}{\sqrt6}K_0^*K_i.
  \]
  Mutual unbiasedness makes its entries unimodular, but does not show that they remain \(24\)th roots of unity. Thus arbitrary unitary normalization does not establish the claimed Butson closure.
- Why it matters: Readers may interpret “a quadruple built from Butson Hadamard matrices” as four Butson bases in an arbitrary ambient coordinate system. The proved theorem only treats a quadruple already represented as
  \(\{I,H_1/\sqrt6,H_2/\sqrt6,H_3/\sqrt6\}\)
  with the three nonidentity matrices in \(\mathrm{BH}(6,24)\).
- Smallest fix: Delete “equivalently,” or define “built from” explicitly to mean the identity-containing form. Do not claim the broader formulation without an additional closure theorem.

### E. The mutual-unbiasedness constant or normalization is mismatched

- Locations: `butson_n24_note.tex:50-58`, `103-109`; `butson_n24_independent.py:49-52`, `122-132`.
- Status: Foreclosed.
- Reason: For unnormalised root-valued columns \(u,w\),
  \[
  \left|\left\langle \frac{u}{\sqrt6},\frac{w}{\sqrt6}\right\rangle\right|^2
  =\frac{|S(w-u)|^2}{36}.
  \]
  Setting this equal to \(1/6\) is exactly the code’s condition \(|S(w-u)|^2=6\). Orthogonality is exactly \(S(w-u)=0\). The identity basis is implicit because every normalized root-valued column has coordinate modulus \(1/\sqrt6\).
- Smallest fix: None.

### F. The transpose indexing error could conceal a matrix-specific completion

- Locations: `butson_n24_note.tex:164-189`; `butson_n24_independent.py:110-116`.
- Status: The per-object correspondence fails, but the universal nonexistence result is foreclosed by transpose closure.
- Reason: A completion of \(H\) omitted at the iteration labelled by \(H\) is tested at the iteration corresponding, up to canonical row sorting, to \(H^T\). The corrected-index recomputation also produced the identical global result.
- Smallest fix: Correct line 114 and regenerate the evidence artefacts so that this extra symmetry argument is unnecessary.

## Final assessment

The narrow theorem stated at `butson_n24_note.tex:73-76` is supported by the exhaustive computation after supplying the elementary transpose-bijection argument. Nevertheless, the current note/checker pair does not satisfy the requested exact line-by-line correspondence because the checker computes \(U(H_1^T)\), and the abstract’s broader “equivalently” claim is not justified. I therefore do not certify the reduction package as clean as written.

"""TASK-111: literature-consistency equivalence-class census on BH(6, N).

For each N in {6, 8, 12, 24}, enumerate dephased row 5-subsets of Z_N-rays
(exactly as `MUBs/butson_exhaust.exhaust_order` does), then quotient these
5-subsets by successively larger dephased-form-preserving groups:

    G0 = S_5^{row}                   (already quotiented by `cliques(...)`)
    G1 = G0 x S_5^{col}              (column permutations on coords 1..5)
    G2 = G1 x <T>                    (matrix transpose)
    G3 = G2 x <K>                    (complex conjugation)

We report the number of orbits of the 5-subset set under G1, G2, G3.
Cross-checks vs. `MUBs/BUTSON_NOTES.md`:

    - N=6:  |Z_6|=340, |M_6|=0     => hadamards computed exactly here.
    - N=8:  432 dephased row-sets, 252 with |U|=0, 180 with |U|=8.
    - N=12: 2184 dephased row-sets, |U|-census {0:1884, 6:240, 12:60}.
    - N=24: 7584 dephased row-sets, |U|-census {0:6852, 6:240, 12:60, 24:360, 40:72}
            (from CLAIM-164, committed at logs/butson_n24_numpy_search.json).

Sub-Butson embeddings BH(6, m) -> BH(6, 24) via v -> (24/m) * v are checked
by counting 5-subsets in N=24 whose entries all lie in (24/m) * Z_m.  For
m | 24 with m < 24, these must equal the BH(6, m) dephased count.

Author family: claude/anthropic (Butson N=24 enumerator authored by google,
independent referee harness by anthropic, independent openai enumerator by
codex; this session performs a structural equivalence-class census).
"""

import argparse
import hashlib
import itertools
import json
import os
import sys
import time
from typing import Dict, Iterable, List, Tuple

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
MUBS_DIR = os.path.join(REPO_ROOT, "MUBs")
SCRIPTS_DIR = os.path.join(REPO_ROOT, "problems", "mub-quad", "scripts")
for p in (MUBS_DIR, SCRIPTS_DIR):
    if p not in sys.path:
        sys.path.insert(0, p)

import butson_exhaust as be  # noqa: E402
from butson_n24_search import ray_classes_numpy  # noqa: E402

DEFAULT_LOG = os.path.join(REPO_ROOT, "problems", "mub-quad", "logs",
                           "butson_bh6_equivalence_classes.json")

D = 6  # matrix order


Ray = Tuple[int, ...]        # (0, v_1, ..., v_5) in (Z/NZ)^6
RowSet = Tuple[Ray, ...]     # sorted tuple of 5 Ray's


def enumerate_dephased_hadamards(n: int) -> Tuple[List[Ray], List[RowSet], Dict[int, int]]:
    """Return (Z-rays, sorted 5-subsets forming BH(6, N), |U|-census).

    Uses the numpy-vectorized ray enumerator for speed on larger N; both
    it and `be.ray_classes` were checked to agree on N in {6, 8, 12} via
    the smoke test in this file's git history.
    """
    if n >= 12:
        zero_rays, mu_rays = ray_classes_numpy(n)
    else:
        zero_rays, mu_rays = be.ray_classes(n)
    # Normalize ordering: sort rays so that downstream be.forward_graph is
    # index-stable across enumerators.
    zero_rays = sorted(zero_rays)
    mu_rays = sorted(mu_rays)
    zero_set = set(zero_rays)
    if not zero_rays:
        return zero_rays, [], {}
    row_graph = be.forward_graph(zero_rays, zero_set, n)
    hadamard_row_sets = be.cliques(row_graph, D - 1)
    mu_set = set(mu_rays)
    u_census: Dict[int, int] = {}
    canon_row_sets: List[RowSet] = []
    zero_row = (0,) * D
    for row_indices in hadamard_row_sets:
        rays = [zero_rays[i] for i in row_indices]
        h_rows = [zero_row] + rays
        u_rays = be.completion_rays(h_rows, mu_rays, mu_set, n)
        u_size = len(u_rays)
        u_census[u_size] = u_census.get(u_size, 0) + 1
        canon_row_sets.append(tuple(sorted(rays)))
    canon_row_sets.sort()
    return zero_rays, canon_row_sets, u_census


def apply_col_perm(rows: RowSet, perm: Tuple[int, ...]) -> RowSet:
    """Permute coordinates 1..5 of each ray by `perm` (a permutation of (1..5))."""
    new_rows = []
    for v in rows:
        new = (0,) + tuple(v[perm[i - 1]] for i in range(1, D))
        new_rows.append(new)
    return tuple(sorted(new_rows))


def transpose_rowset(rows: RowSet) -> RowSet:
    """Return the row 5-subset of the transposed dephased matrix.

    H has 6 rows: row_0 = (0)^6, and rows_{1..5} = rows[0..4] (in sorted order).
    H^T has row_j = column_j of H = (0, rows[0][j], rows[1][j], ..., rows[4][j]).
    Only j = 1..5 give nonzero rows; row 0 of H^T is (0)^6.
    """
    assert len(rows) == D - 1
    new_rows = []
    for j in range(1, D):
        col_j = (0,) + tuple(rows[i][j] for i in range(D - 1))
        new_rows.append(col_j)
    return tuple(sorted(new_rows))


def conjugate_rowset(rows: RowSet, n: int) -> RowSet:
    """Complex-conjugation: replace zeta^v by zeta^{-v}, i.e., v -> -v mod N."""
    new_rows = []
    for v in rows:
        new = (0,) + tuple((-v[i]) % n for i in range(1, D))
        new_rows.append(new)
    return tuple(sorted(new_rows))


def orbit_canonical(rows: RowSet, n: int, group: str) -> RowSet:
    """Canonical representative under the chosen group.

    group ∈ {"col", "col_T", "col_T_K"}.
    """
    include_T = "T" in group
    include_K = "K" in group
    perms = list(itertools.permutations(range(1, D)))  # S_5 on coords 1..5
    best: RowSet | None = None
    seeds: List[RowSet] = [rows]
    if include_T:
        seeds.append(transpose_rowset(rows))
    if include_K:
        conj_seeds = [conjugate_rowset(s, n) for s in seeds]
        seeds.extend(conj_seeds)
    for seed in seeds:
        for perm in perms:
            candidate = apply_col_perm(seed, perm)
            if best is None or candidate < best:
                best = candidate
    assert best is not None
    return best


def orbit_class_count(rowsets: List[RowSet], n: int, group: str) -> int:
    seen = set()
    for rs in rowsets:
        seen.add(orbit_canonical(rs, n, group))
    return len(seen)


def sub_butson_embedding_count(rowsets: List[RowSet], n: int, m: int) -> int:
    """Count row-sets in BH(6, N) whose entries all lie in the sub-lattice (N/m) * Z/mZ.

    These correspond to BH(6, m) matrices under the embedding zeta_m -> zeta_N^{N/m}.
    Precondition: m divides n.
    """
    assert n % m == 0
    scale = n // m
    hits = 0
    for rs in rowsets:
        if all(entry % scale == 0 for row in rs for entry in row):
            hits += 1
    return hits


def compute_bh6_census(orders: List[int], with_transpose_conjugation: bool,
                       progress: bool = False) -> Dict:
    """Compute dephased-count and equivalence-class census for each N in `orders`."""
    results: Dict[str, Dict] = {}
    order_metadata: Dict[str, Dict] = {}
    for n in orders:
        if progress:
            print(f"  computing BH(6, {n}) ...", flush=True)
        started = time.time()
        _, rowsets, u_census = enumerate_dephased_hadamards(n)
        enum_seconds = time.time() - started
        entry: Dict = {
            "dephased_row_sets": len(rowsets),
            "u_census": {str(k): u_census[k] for k in sorted(u_census)},
        }
        if rowsets:
            entry["equivalence_classes"] = {}
            for group_label in ("col", "col_T", "col_T_K"):
                if group_label != "col" and not with_transpose_conjugation:
                    continue
                t0 = time.time()
                count = orbit_class_count(rowsets, n, group_label)
                entry["equivalence_classes"][group_label] = count
                entry.setdefault("timing_seconds", {})[f"orbit_{group_label}"] = round(
                    time.time() - t0, 3
                )
        order_metadata[str(n)] = {"enumeration_seconds": round(enum_seconds, 3)}
        results[str(n)] = entry
        if progress:
            print(f"    dephased={entry['dephased_row_sets']}", flush=True)
            if "equivalence_classes" in entry:
                print(f"    equiv classes: {entry['equivalence_classes']}", flush=True)
    return {
        "orders": orders,
        "with_transpose_conjugation": with_transpose_conjugation,
        "per_order": results,
        "timing": order_metadata,
    }


def sub_butson_check(orders_for_embedding: List[int], n_top: int,
                     rowsets_top: List[RowSet],
                     per_order_dephased: Dict[int, int]) -> Dict[str, int]:
    """Sanity: for each m | n_top with m < n_top, count top-N sets that come from BH(6,m)."""
    out: Dict[str, int] = {}
    for m in orders_for_embedding:
        if m >= n_top or n_top % m != 0:
            continue
        count = sub_butson_embedding_count(rowsets_top, n_top, m)
        out[str(m)] = count
        expected = per_order_dephased.get(m)
        if expected is not None:
            out[f"{m}_matches_bh6_{m}"] = (count == expected)
    return out


def compute_full_summary(orders: List[int], top_order: int,
                         progress: bool = False) -> Dict:
    """Full pipeline: census + embedding cross-check."""
    census = compute_bh6_census(orders, with_transpose_conjugation=True,
                                progress=progress)

    # Redo enumeration only for the top order to keep the row-set list for
    # embedding cross-check.  This is cheap; keeps compute_bh6_census pure.
    if progress:
        print(f"  cross-check: embedding lower orders into BH(6, {top_order})",
              flush=True)
    _, rowsets_top, _ = enumerate_dephased_hadamards(top_order)
    per_order_dephased: Dict[int, int] = {
        int(k): v["dephased_row_sets"] for k, v in census["per_order"].items()
    }
    embed = sub_butson_check(orders, top_order, rowsets_top, per_order_dephased)
    census["sub_butson_embeddings_into_N"] = {
        "top_order": top_order,
        "embeddings": embed,
    }
    return census


def deterministic_content(payload: Dict) -> Dict:
    """Return a copy with wall-clock fields stripped for deterministic hashing."""
    stripped = json.loads(json.dumps(payload, sort_keys=True))
    stripped.pop("wall_seconds", None)
    stripped.pop("timing", None)
    stripped.pop("deterministic_content_sha256", None)
    for entry in stripped.get("per_order", {}).values():
        entry.pop("timing_seconds", None)
    return stripped


def content_sha256(payload: Dict) -> str:
    stripped = deterministic_content(payload)
    canonical = json.dumps(stripped, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def write_log(payload: Dict, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def cmd_compute(args) -> None:
    started = time.time()
    orders = [int(x) for x in args.orders.split(",")]
    orders = sorted(set(orders))
    top_order = args.top_order
    summary = compute_full_summary(orders, top_order, progress=args.progress)
    summary["wall_seconds"] = round(time.time() - started, 3)
    summary["deterministic_content_sha256"] = content_sha256(summary)
    write_log(summary, args.log)
    print(f"Wrote {args.log}", flush=True)
    print(f"deterministic_content_sha256 = {summary['deterministic_content_sha256']}",
          flush=True)


def cmd_verify(args) -> None:
    with open(args.log) as f:
        recorded = json.load(f)
    expected_hash = recorded.get("deterministic_content_sha256")
    computed = content_sha256(recorded)
    assert expected_hash == computed, (
        f"deterministic_content_sha256 mismatch: recorded={expected_hash}, "
        f"computed={computed}"
    )
    # Regression: cross-check counts against BUTSON_NOTES.md and CLAIM-164.
    per_order = recorded["per_order"]
    assertions = {
        "6":  {"dephased_row_sets": None},   # not published in BUTSON_NOTES; compute-only
        "8":  {"dephased_row_sets": 432},
        "12": {"dephased_row_sets": 2184},
        "24": {"dephased_row_sets": 7584},
    }
    u_expected = {
        "8":  {"0": 252, "8": 180},
        "12": {"0": 1884, "6": 240, "12": 60},
        "24": {"0": 6852, "6": 240, "12": 60, "24": 360, "40": 72},
    }
    for n, exp in assertions.items():
        got = per_order.get(n)
        assert got is not None, f"missing per_order[{n}]"
        if exp["dephased_row_sets"] is not None:
            assert got["dephased_row_sets"] == exp["dephased_row_sets"], (
                f"BH(6, {n}) count mismatch: got {got['dephased_row_sets']}, "
                f"expected {exp['dephased_row_sets']}"
            )
        if n in u_expected:
            assert got["u_census"] == u_expected[n], (
                f"BH(6, {n}) |U|-census mismatch: got {got['u_census']}, "
                f"expected {u_expected[n]}"
            )
    # Sub-Butson embeddings into N=24: BH(6, 8) embed count = 432, BH(6, 12) = 2184.
    embed = recorded["sub_butson_embeddings_into_N"]["embeddings"]
    if "8" in embed:
        assert embed["8"] == 432, f"BH(6, 8) embedding count: {embed['8']}"
        assert embed["8_matches_bh6_8"] is True
    if "12" in embed:
        assert embed["12"] == 2184, f"BH(6, 12) embedding count: {embed['12']}"
        assert embed["12_matches_bh6_12"] is True
    print("VERIFIED:", args.log, "sha256 =", expected_hash)


def cmd_selftest(args) -> None:
    """Micro-test: apply group actions on a hand-designed row-set for BH(6, 8)."""
    # For BH(6, 8), pick a canonical row-set: identity 6x6 diagonal Fourier
    # would need entries in {0,1,2,3,4,5,6,7}.  Instead we just check that
    # applying the identity permutation leaves the row-set fixed and that
    # transpose is idempotent up to double-application when the matrix is
    # symmetric.  Concretely, pick the first BH(6, 8) row-set.
    _, rowsets, _ = enumerate_dephased_hadamards(8)
    assert rowsets, "BH(6, 8) should be non-empty"
    rs = rowsets[0]
    # Identity permutation is idempotent.
    perm_id = tuple(range(1, D))
    assert apply_col_perm(rs, perm_id) == rs
    # Transpose applied twice on the underlying matrix is identity.  Check
    # that transposing twice yields the original (as a sorted 5-tuple).
    tt = transpose_rowset(transpose_rowset(rs))
    assert tt == rs, f"transpose ** 2 != id: {tt} vs {rs}"
    # Conjugation applied twice is identity.
    kk = conjugate_rowset(conjugate_rowset(rs, 8), 8)
    assert kk == rs, f"conj ** 2 != id: {kk} vs {rs}"
    print("SELF-TEST PASSED: identity, transpose^2, conj^2 fix the first BH(6, 8) row-set.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_compute = sub.add_parser("compute")
    p_compute.add_argument("--orders", default="6,8,12,24",
                           help="Comma-separated list of N values")
    p_compute.add_argument("--top-order", type=int, default=24,
                           help="Order used for sub-Butson embedding cross-check")
    p_compute.add_argument("--log", default=DEFAULT_LOG)
    p_compute.add_argument("--progress", action="store_true")
    p_compute.set_defaults(func=cmd_compute)

    p_verify = sub.add_parser("verify")
    p_verify.add_argument("--log", default=DEFAULT_LOG)
    p_verify.set_defaults(func=cmd_verify)

    p_self = sub.add_parser("self-test")
    p_self.set_defaults(func=cmd_selftest)

    args = parser.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())

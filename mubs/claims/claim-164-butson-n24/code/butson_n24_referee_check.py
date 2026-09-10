#!/usr/bin/env python3
"""Referee cross-checks for CLAIM-164 (Butson N=24 MUB quadruple exhaustion).

Independent verification harness written by the anthropic referee session of
2026-08-21. It re-derives the CLAIM-164 census through code paths that are
independent of `butson_n24_search.py` wherever feasible, and re-verifies the
load-bearing exclusion (zero mutually unbiased completion-basis pairs) with
exact arithmetic in Z[zeta_24] via `MUBs/mub_common.Cyclo`, a representation
(coefficients mod x^N - 1, reduced by Phi_N on demand) disjoint from the
search's power-basis reduction tensor.

Stages (run from the repo root; state is passed through /tmp):

  --stage 1  Foundations.
             (A) Hand-coded Phi_24 = x^8 - x^4 + 1 reduction table compared
                 entry-by-entry with `cyclotomic_reductions(24)` (sympy), and
                 numerically against zeta_24^k in complex128.
             (B) Set equality (not just counts) of `ray_classes_numpy(n)`
                 against the exact pure-Python `butson_exhaust.ray_classes(n)`
                 at n = 8 and n = 12.
             (D) Full-pipeline regression `butson_exhaust.exhaust_order(12)`
                 against the independently authored EXPECTED census.
  --stage 2  N = 24 rebuild.
             (C) Float classification of all 24^5 dephased tuples in
                 complex128 with margin analysis; near-boundary cases (if any)
                 adjudicated exactly with Cyclo; sets compared to
                 `ray_classes_numpy(24)`.
             Row graph built with independent numpy key arithmetic (not
             `butson_exhaust.forward_graph`), 5-cliques via
             `butson_exhaust.cliques`; count compared to the committed 7584.
             U(H1) computed for every H1 with an independent vectorised
             membership test; |U| census compared to the committed log;
             `butson_exhaust.completion_rays` cross-checked on samples.
  --stage 3  Bases and exact exclusion.
             Completion bases (6-cliques) for every H1 with U != 0; the full
             |U|/bases census compared to the committed log. Every base is
             checked exactly (Cyclo): orthonormality and unbiasedness to its
             H1. Every candidate pair of bases (H2, H3) sharing an H1 is
             exactly refuted: a cross column pair with |<x, y>|^2 != 6 in
             Z[zeta_24] is exhibited. Writes the combined deterministic log.

Exit status is nonzero on any mismatch. The combined log is written to
`problems/mub-quad/logs/butson_n24_referee_check.json`; wall-clock timings
are excluded from the deterministic-content SHA-256.
"""

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "MUBs")
sys.path.insert(0, "problems/mub-quad/scripts")

import mub_common as mc
import butson_exhaust as be
from butson_n24_vanishing_sum import cyclotomic_reductions
from butson_n24_search import ray_classes_numpy

N = 24
D = 6
STATE = "/tmp/butson_n24_referee_state.npz"
PARTIAL = "/tmp/butson_n24_referee_partial.json"
LOG = "problems/mub-quad/logs/butson_n24_referee_check.json"

COMMITTED = {
    "zero_rays": 8350,
    "mu_rays": 19320,
    "hadamards": 7584,
    "total_completions": 1920,
    "max_u": 40,
    "max_completions": 4,
    "types": {"0/0": 6852, "12/4": 60, "24/4": 360, "40/0": 72, "6/1": 240},
}

POWS = (N ** np.arange(D)).astype(np.int64)


def keys_of(arr):
    """Injective int64 key of exponent vectors (last axis length 6)."""
    return arr.astype(np.int64) @ POWS


def hand_reductions_24():
    """reps[k] for x^k mod Phi_24, Phi_24 = x^8 - x^4 + 1 hard-coded."""
    phi = [1, 0, 0, 0, -1, 0, 0, 0, 1]  # low degree first
    reps = []
    cur = [0] * 8
    cur[0] = 1
    for _ in range(N):
        reps.append(cur[:])
        new = [0] * 8
        for i in range(7):
            new[i + 1] = cur[i]
        top = cur[7]
        if top:
            for i in range(8):
                new[i] -= top * phi[i]
        cur = new
    return reps


def stage1():
    out = {}
    # (A) reduction table: hand-coded vs committed generator vs numeric.
    mine = hand_reductions_24()
    theirs = cyclotomic_reductions(24)
    assert [list(r) for r in theirs] == mine, "reduction table mismatch"
    assert mine[0] == [1, 0, 0, 0, 0, 0, 0, 0], "basis must start at 1"
    zeta = np.exp(2j * np.pi / N)
    basis = zeta ** np.arange(8)
    err = max(
        abs((np.array(mine[k]) * basis).sum() - zeta ** k) for k in range(N)
    )
    assert err < 1e-12, f"numeric reduction error {err}"
    assert mc.cyclotomic_poly(24) == [1, 0, 0, 0, -1, 0, 0, 0, 1], \
        "mub_common Phi_24 mismatch"
    out["A_reduction_table"] = {"agrees": True, "numeric_err": float(err)}

    # (B) ray classes: numpy path vs exact Cyclo path, full set equality.
    for n in (8, 12):
        z_np, m_np = ray_classes_numpy(n)
        z_ex, m_ex = be.ray_classes(n)
        assert set(z_np) == set(z_ex), f"Z set mismatch at n={n}"
        assert set(m_np) == set(m_ex), f"M set mismatch at n={n}"
        out[f"B_ray_sets_n{n}"] = {"Z": len(z_ex), "M": len(m_ex),
                                   "sets_equal": True}

    # (D) full pipeline regression at N=12 (raises on census mismatch).
    r12 = be.exhaust_order(12)
    assert r12["verdict"] == "no_quadruple"
    out["D_pipeline_n12"] = {k: r12[k] for k in
                            ("hadamards", "total_completions", "types",
                             "verdict")}
    return out


def stage2():
    out = {}
    t0 = time.time()
    Z, M = ray_classes_numpy(N)
    assert len(Z) == COMMITTED["zero_rays"] and len(M) == COMMITTED["mu_rays"]
    Zarr = np.array(Z, dtype=np.int16)
    Marr = np.array(M, dtype=np.int16)
    zkeys = np.sort(keys_of(Zarr))
    mkeys = np.sort(keys_of(Marr))

    # (C) independent float classification of all 24^5 dephased tuples.
    zeta_pow = np.exp(2j * np.pi * np.arange(N) / N)
    inner = np.stack(
        np.meshgrid(*([np.arange(N)] * 4), indexing="ij"), axis=-1
    ).reshape(-1, 4).astype(np.int16)
    cnt = inner.shape[0]
    z_float = 0
    m_float = 0
    margin_z = np.inf   # min |S| among non-zero-sum tuples
    margin_m = np.inf   # min ||S|^2 - 6| among non-mu tuples
    suspects = 0        # tuples inside the ambiguous float band
    for v1 in range(N):
        vec = np.concatenate(
            [np.zeros((cnt, 1), np.int16),
             np.full((cnt, 1), v1, np.int16), inner], axis=1)
        S = zeta_pow[vec].sum(axis=1)
        a = np.abs(S)
        q = np.abs(a * a - 6.0)
        zm = a < 1e-6
        mm = q < 1e-6
        suspects += int(((a > 1e-9) & (a < 1e-3)).sum())
        suspects += int(((q > 1e-9) & (q < 1e-3)).sum())
        z_float += int(zm.sum())
        m_float += int(mm.sum())
        if (~zm).any():
            margin_z = min(margin_z, float(a[~zm].min()))
        if (~mm).any():
            margin_m = min(margin_m, float(q[~mm].min()))
        # membership must agree with the exact sets tuple-for-tuple
        k = keys_of(vec)
        in_z = zkeys[np.searchsorted(zkeys, k).clip(0, len(zkeys) - 1)] == k
        in_m = mkeys[np.searchsorted(mkeys, k).clip(0, len(mkeys) - 1)] == k
        assert (zm == in_z).all(), f"Z float/exact split at v1={v1}"
        assert (mm == in_m).all(), f"M float/exact split at v1={v1}"
    assert suspects == 0, f"{suspects} tuples in ambiguous float band"
    assert z_float == len(Z) and m_float == len(M)
    out["C_float_recount"] = {
        "Z": z_float, "M": m_float, "sets_equal": True,
        "min_margin_nonzero_absS": margin_z,
        "min_margin_nonmu_abs_q": margin_m,
        "ambiguous_band_tuples": suspects,
    }

    # Row graph with independent numpy arithmetic; cliques via be.cliques.
    nz = len(Z)
    adjacency = [0] * nz
    edges = 0
    chunk = 512
    for lo in range(0, nz, chunk):
        hi = min(lo + chunk, nz)
        diffs = (Zarr[None, :, :] - Zarr[lo:hi, None, :]) % N  # (c, nz, 6)
        k = keys_of(diffs)
        pos = np.searchsorted(zkeys, k).clip(0, len(zkeys) - 1)
        member = zkeys[pos] == k
        for i in range(lo, hi):
            row = member[i - lo]
            js = np.nonzero(row)[0]
            js = js[js > i]
            bits = 0
            for j in js:
                bits |= 1 << int(j)
            adjacency[i] = bits
            edges += len(js)
    h1_sets = be.cliques(adjacency, 5)
    assert len(h1_sets) == COMMITTED["hadamards"], \
        f"hadamards {len(h1_sets)} != {COMMITTED['hadamards']}"
    out["row_graph"] = {"edges": edges, "hadamards": len(h1_sets)}

    # Exact Cyclo row-orthogonality on a deterministic sample of H1s.
    zero_row = (0,) * D
    sampled = 0
    for idx in range(0, len(h1_sets), 37):
        rows = [zero_row] + [Z[i] for i in h1_sets[idx]]
        h = mc.root_matrix(N, rows)
        hh = mc.cyclo_mat_mul_dagger_left(h, h)
        for a in range(D):
            for b in range(D):
                tgt = mc.Cyclo.integer(N, D if a == b else 0)
                assert mc.cyclo_equal(hh[a][b], tgt), \
                    f"H1 #{idx} rows not exactly orthogonal"
        sampled += 1
    out["exact_h1_sample"] = {"stride": 37, "checked": sampled, "ok": True}

    # U(H1) for every H1, via independent vectorised membership test.
    u_census = {}
    u_lists = []
    cross_checked = set()
    for idx, row_set in enumerate(h1_sets):
        rows = [zero_row] + [Z[i] for i in row_set]
        cols = np.array(list(zip(*rows)), dtype=np.int16)  # (6, 6)
        diffs = (Marr[None, :, :] - cols[1:, None, :]) % N  # skip all-0 col
        k = keys_of(diffs)
        pos = np.searchsorted(mkeys, k).clip(0, len(mkeys) - 1)
        member = (mkeys[pos] == k).all(axis=0)
        u_idx = np.nonzero(member)[0]
        u = len(u_idx)
        u_census[u] = u_census.get(u, 0) + 1
        u_lists.append([int(x) for x in u_idx])
        if u and u not in cross_checked:
            mine = [M[i] for i in u_idx]
            ref = be.completion_rays(rows, M, set(M), N)
            assert mine == ref, f"U(H1) mismatch vs completion_rays at #{idx}"
            cross_checked.add(u)
    expected_u = {}
    for key, count in COMMITTED["types"].items():
        u = int(key.split("/")[0])
        expected_u[u] = expected_u.get(u, 0) + count
    assert u_census == expected_u, f"|U| census {u_census} != {expected_u}"
    out["u_census"] = {str(k): v for k, v in sorted(u_census.items())}
    out["u_cross_checked_classes"] = sorted(cross_checked)
    out["seconds"] = round(time.time() - t0, 3)

    np.savez_compressed(
        STATE,
        Z=Zarr, M=Marr,
        h1=np.array(h1_sets, dtype=np.int32),
        u_flat=np.array([i for lst in u_lists for i in lst], dtype=np.int32),
        u_len=np.array([len(lst) for lst in u_lists], dtype=np.int32),
    )
    return out


def stage3():
    out = {}
    t0 = time.time()
    st = np.load(STATE)
    Z = [tuple(int(x) for x in r) for r in st["Z"]]
    M = [tuple(int(x) for x in r) for r in st["M"]]
    h1_sets = [tuple(int(x) for x in r) for r in st["h1"]]
    u_len = st["u_len"]
    u_flat = st["u_flat"]
    offsets = np.concatenate([[0], np.cumsum(u_len)])
    zero_set = set(Z)
    mu_set = set(M)
    zero_row = (0,) * D
    six = mc.Cyclo.integer(N, D)

    census = {}
    total_bases = 0
    max_completions = 0
    n_nonzero_u = 0
    n_with_base = 0
    pairs_refuted = 0
    exact_base_checks = 0
    for idx, row_set in enumerate(h1_sets):
        u = int(u_len[idx])
        if u == 0:
            census["0/0"] = census.get("0/0", 0) + 1
            continue
        n_nonzero_u += 1
        rays = [M[i] for i in u_flat[offsets[idx]:offsets[idx + 1]]]
        graph = be.forward_graph(rays, zero_set, N)
        bases = be.cliques(graph, D)
        key = "%d/%d" % (u, len(bases))
        census[key] = census.get(key, 0) + 1
        total_bases += len(bases)
        max_completions = max(max_completions, len(bases))
        if not bases:
            continue
        n_with_base += 1
        rows = [zero_row] + [Z[i] for i in row_set]
        h1 = mc.root_matrix(N, rows)
        for base in bases:
            exp = [[rays[j][i] for j in base] for i in range(D)]
            b = mc.root_matrix(N, exp)
            assert mc.exact_orthonormal(b, D), \
                f"base not exactly orthonormal at H1 #{idx}"
            assert mc.exact_mu_pair(h1, b, D), \
                f"base not exactly unbiased to H1 #{idx}"
            exact_base_checks += 1
        for a in range(len(bases)):
            for c in range(a + 1, len(bases)):
                assert not be.compatible_completions(
                    bases[a], bases[c], rays, mu_set, N)
                witness = None
                for i in bases[a]:
                    for j in bases[c]:
                        val = mc.Cyclo(N)
                        for s in range(D):
                            val = val + mc.Cyclo.root(
                                N, rays[j][s] - rays[i][s])
                        if not mc.cyclo_equal(val * val.conj(), six):
                            witness = (i, j)
                            break
                    if witness:
                        break
                assert witness is not None, \
                    f"pair at H1 #{idx} has all 36 cross pairs unbiased!"
                pairs_refuted += 1

    assert census == COMMITTED["types"], \
        f"census {census} != committed {COMMITTED['types']}"
    assert total_bases == COMMITTED["total_completions"]
    assert max_completions == COMMITTED["max_completions"]
    out["bases_census"] = census
    out["total_completions"] = total_bases
    out["hadamards_with_nonzero_U"] = n_nonzero_u
    out["hadamards_with_at_least_one_base"] = n_with_base
    out["exact_base_checks"] = exact_base_checks
    out["candidate_pairs_exactly_refuted"] = pairs_refuted
    out["seconds"] = round(time.time() - t0, 3)
    return out


def strip_timings(obj):
    if isinstance(obj, dict):
        return {k: strip_timings(v) for k, v in obj.items()
                if k != "seconds"}
    if isinstance(obj, list):
        return [strip_timings(v) for v in obj]
    return obj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, required=True, choices=(1, 2, 3))
    args = ap.parse_args()

    partial = {}
    if os.path.exists(PARTIAL):
        with open(PARTIAL) as fh:
            partial = json.load(fh)

    result = {1: stage1, 2: stage2, 3: stage3}[args.stage]()
    partial[f"stage{args.stage}"] = result
    with open(PARTIAL, "w") as fh:
        json.dump(partial, fh, indent=2, sort_keys=True)
    print(f"stage {args.stage} OK")

    if args.stage == 3:
        partial["committed_reference"] = COMMITTED
        payload = json.dumps(strip_timings(partial), sort_keys=True,
                             separators=(",", ":")).encode()
        partial["deterministic_sha256"] = hashlib.sha256(payload).hexdigest()
        with open(LOG, "w") as fh:
            json.dump(partial, fh, indent=2, sort_keys=True)
        print(f"wrote {LOG}")
        print(f"deterministic SHA-256: {partial['deterministic_sha256']}")


if __name__ == "__main__":
    main()

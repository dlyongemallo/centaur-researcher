#!/usr/bin/env python3
"""TASK-84 pilot: PB-SAT enumeration of Butson vanishing-sum tuples.

For each cyclotomic order N, count the dephased root-valued 6-tuples

    Z_N = { (0, v_1, ..., v_5) in (Z/NZ)^6 :  sum_j zeta_N^{v_j} = 0 }
    M_N = { (0, v_1, ..., v_5) in (Z/NZ)^6 : |sum_j zeta_N^{v_j}|^2 = 6 }

directly used by the Butson $N$-th-root MUB exhaustion of
`MUBs/BUTSON_NOTES.md`. The vanishing condition is linear in the one-hot
slot indicators, so it reduces to $\varphi(N)$ pseudo-boolean equalities on
the cyclotomic basis of Z[zeta_N]; the pilot encodes each equality via a
CardEnc.equals cardinality constraint with literal flips (all reduction
coefficients lie in {-1, 0, 1} for N in {6, 8, 12, 24}).

The pilot cross-checks the SAT enumerator against a self-contained brute-
force enumerator that evaluates each tuple exactly in the cyclotomic basis
(and thus matches `MUBs/butson_exhaust.ray_classes` up to the CardEnc/basis
choice). |M_N| is quadratic in the one-hot indicators; the pilot enumerates
it by brute force only.

Usage:
    python3 butson_n24_vanishing_sum.py --verify
    python3 butson_n24_vanishing_sum.py --brute N        # counts Z_N, M_N
    python3 butson_n24_vanishing_sum.py --sat-z N        # SAT enumerates Z_N
    python3 butson_n24_vanishing_sum.py --self-test

CLAIM tier: T1 pilot; SAT enumerator and brute-force enumerator both cold-
verified against the BUTSON_NOTES.md census at N in {6, 8, 12}. N = 24 is
recorded with an explicit budget and outcome; a heavier N = 24 pass belongs
in `run_long`.
"""

import argparse
import hashlib
import itertools
import json
import os
import sys
import time
from fractions import Fraction

# Reference counts transcribed verbatim from MUBs/BUTSON_NOTES.md.
# The two hand-copied dictionaries below are the ONLY data the verifier
# trusts as ground truth; the SAT and brute-force enumerators must agree
# with them at every listed N.
REF_Z = {
    1: 0, 2: 10, 3: 30, 4: 100, 5: 0, 6: 340, 7: 0, 8: 640,
    9: 270, 10: 1090, 11: 0, 12: 1930,
}
REF_M = {
    1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1680,
    9: 0, 10: 0, 11: 0, 12: 3600,
}

DIM = 6

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DEFAULT_LOG = os.path.join(REPO_ROOT, "problems", "mub-quad", "logs",
                           "butson_n24_vanishing_sum.json")


def cyclotomic_reductions(n):
    """Return reps[k] = coeff vector of x^k mod Phi_n(x) for k = 0..n-1.

    Each reps[k] is a length-phi(n) list of integers so that
    x^k = sum_i reps[k][i] * x^i (mod Phi_n(x)).
    """
    from sympy import cyclotomic_poly, Poly, symbols
    x = symbols("x")
    phi_poly = Poly(cyclotomic_poly(n, x), x)
    phi_n = phi_poly.degree()
    coeffs = list(reversed(phi_poly.all_coeffs()))
    assert coeffs[-1] == 1, "cyclotomic leading coeff must be 1"
    reps = []
    cur = [0] * phi_n
    cur[0] = 1
    for _ in range(n):
        reps.append(cur[:])
        new = [0] * phi_n
        for i in range(phi_n - 1):
            new[i + 1] = cur[i]
        top = cur[phi_n - 1]
        if top:
            for i in range(phi_n):
                new[i] -= top * coeffs[i]
        cur = new
    return reps


def brute_counts(n):
    """Enumerate |Z_N| and |M_N| directly in the cyclotomic basis (pure Python)."""
    reps = cyclotomic_reductions(n)
    phi_n = len(reps[0])
    zc = 0
    mc = 0
    for tail in itertools.product(range(n), repeat=DIM - 1):
        vec = (0,) + tail
        total = [0] * phi_n
        for k in vec:
            for i in range(phi_n):
                total[i] += reps[k][i]
        if not any(total):
            zc += 1
        norm = [0] * phi_n
        for a in vec:
            for b in vec:
                d = (a - b) % n
                for i in range(phi_n):
                    norm[i] += reps[d][i]
        if norm[0] == DIM and not any(norm[1:]):
            mc += 1
    return zc, mc


def brute_counts_numpy(n, batch_outer=None):
    """Vectorized brute-force |Z_N|, |M_N|. Batches by the outermost tail
    coordinate to keep memory bounded; independent numpy code path from
    brute_counts() serving as its own cross-check."""
    import numpy as np
    reps_list = cyclotomic_reductions(n)
    reps = np.array(reps_list, dtype=np.int8)  # shape (n, phi_n)
    phi_n = reps.shape[1]
    slots_tail = DIM - 1  # v_1..v_5
    outer_choices = list(range(n))
    if batch_outer is None:
        batch_outer = 1  # one v_1 value at a time
    zc = 0
    mc = 0
    # tail = (v_1, v_2, v_3, v_4, v_5). We iterate over v_1 and vectorize
    # v_2..v_5. That is n^4 = 24^4 = 331,776 per batch at N=24; memory
    # for the 6x6 pairwise-difference tensor is 331_776 * 36 * phi_n bytes
    # ~= 76 MB for phi_n=8, easily inside 4 GB.
    inner_shape = (n,) * (slots_tail - 1)
    inner_grid = np.stack(np.meshgrid(*([np.arange(n)] * (slots_tail - 1)),
                                      indexing="ij"), axis=-1)
    inner_flat = inner_grid.reshape(-1, slots_tail - 1).astype(np.int16)
    inner_count = inner_flat.shape[0]

    for v1 in outer_choices:
        # vec = (0, v1, inner[:, 0], inner[:, 1], inner[:, 2], inner[:, 3])
        v1_col = np.full((inner_count, 1), v1, dtype=np.int16)
        z0_col = np.zeros((inner_count, 1), dtype=np.int16)
        vec = np.concatenate([z0_col, v1_col, inner_flat], axis=1)  # (M, 6)

        # Total sum in cyclotomic basis: sum reps[v_j]
        # reps: (n, phi_n), vec: (M, 6) -> gather reps[vec] (M, 6, phi_n)
        gathered = reps[vec]  # (M, 6, phi_n)
        totals = gathered.sum(axis=1)  # (M, phi_n)
        zc += int(np.all(totals == 0, axis=1).sum())

        # |S|^2 = sum_{a,b} reps[(v_a - v_b) mod n]
        diffs = (vec[:, :, None] - vec[:, None, :]) % n  # (M, 6, 6)
        # Gather reps[diffs] -> (M, 6, 6, phi_n) then sum over (1,2)
        diffs_reps = reps[diffs]  # (M, 6, 6, phi_n)
        norms = diffs_reps.sum(axis=(1, 2))  # (M, phi_n)
        m_hit = (norms[:, 0] == DIM) & (norms[:, 1:].sum(axis=1) == 0) & \
                np.all(norms[:, 1:] == 0, axis=1)
        mc += int(m_hit.sum())

    return zc, mc


def sat_enumerate_Z(n, cap=None):
    """SAT-enumerate Z_N via CardEnc.equals; return (count, capped_flag, seconds)."""
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Solver
    from pysat.formula import IDPool

    reps = cyclotomic_reductions(n)
    phi_n = len(reps[0])
    slots = DIM - 1  # slot 0 is fixed at v_0 = 0

    pool = IDPool()
    # x[j][k] for j in 0..slots-1 (representing v_{j+1}), k in 0..n-1
    x = [[pool.id(f"x_{j}_{k}") for k in range(n)] for j in range(slots)]

    solver = Solver(name="cadical195")

    # Exactly-one per slot: reps[0] = (1, 0, ..., 0), so the sum
    # over one-hot k of reps[k][i] gives the coeff of x^i in that slot.
    for j in range(slots):
        eq = CardEnc.equals(lits=x[j], bound=1, top_id=pool.top,
                            encoding=EncType.pairwise)
        pool.top = max(pool.top, eq.nv)
        for cl in eq.clauses:
            solver.add_clause(cl)

    # PB equality per basis dim: sum_{j,k} reps[k][i] * x[j][k] = -[i==0]
    # (the -1 comes from the fixed slot 0 contributing +1 to dim 0).
    for i in range(phi_n):
        target_rhs = -1 if i == 0 else 0
        pos_lits = []
        neg_lits = []
        for j in range(slots):
            for k in range(n):
                c = reps[k][i]
                if c == 1:
                    pos_lits.append(x[j][k])
                elif c == -1:
                    neg_lits.append(x[j][k])
        # Flip transform: sum_P x + sum_N (1 - x) = target_rhs + |N|
        lits = list(pos_lits) + [-v for v in neg_lits]
        bound = target_rhs + len(neg_lits)
        if bound < 0 or bound > len(lits):
            solver.delete()
            return 0, False, 0.0  # infeasible on this dim
        eq = CardEnc.equals(lits=lits, bound=bound, top_id=pool.top,
                            encoding=EncType.seqcounter)
        pool.top = max(pool.top, eq.nv)
        for cl in eq.clauses:
            solver.add_clause(cl)

    # AllSAT loop: enumerate satisfying assignments over just the one-hot vars.
    count = 0
    capped = False
    t0 = time.perf_counter()
    while solver.solve():
        model = solver.get_model()
        # Extract picked v_1..v_5
        picks = []
        picked_lits = []
        for j in range(slots):
            picked = None
            for k in range(n):
                if x[j][k] in model:
                    assert picked is None, "slot not one-hot"
                    picked = k
                    picked_lits.append(x[j][k])
            assert picked is not None
            picks.append(picked)
        count += 1
        # Block this exact one-hot pattern
        solver.add_clause([-l for l in picked_lits])
        if cap is not None and count >= cap:
            capped = True
            break
    dt = time.perf_counter() - t0
    solver.delete()
    return count, capped, dt


def self_test():
    """Cross-check SAT, pure-Python brute, and numpy brute at N in {6, 8}."""
    ok = True
    for n in (6, 8):
        zb, mb = brute_counts(n)
        zn, mn = brute_counts_numpy(n)
        zs, capped, _ = sat_enumerate_Z(n)
        assert zb == REF_Z[n], f"brute Z_{n}={zb} != ref {REF_Z[n]}"
        assert mb == REF_M[n], f"brute M_{n}={mb} != ref {REF_M[n]}"
        assert zn == zb and mn == mb, f"numpy vs pure disagreement at N={n}"
        assert zs == REF_Z[n], f"SAT Z_{n}={zs} != ref {REF_Z[n]}"
        assert not capped
        print(f"  N={n}: brute Z={zb} M={mb}; numpy Z={zn} M={mn}; SAT Z={zs}  -- ok")
    # Phi_24 reduction integrity: x^12 = -1
    reps24 = cyclotomic_reductions(24)
    assert reps24[12] == [-1, 0, 0, 0, 0, 0, 0, 0], f"x^12 reduction wrong: {reps24[12]}"
    assert reps24[0] == [1, 0, 0, 0, 0, 0, 0, 0]
    assert reps24[8] == [-1, 0, 0, 0, 1, 0, 0, 0], f"x^8 reduction wrong: {reps24[8]}"
    for k in range(24):
        for c in reps24[k]:
            assert c in (-1, 0, 1), f"reps24[{k}] has coeff {c} outside +/-1"
    print("  Phi_24 reduction: x^0, x^8, x^12 correct; all coefficients in {-1, 0, 1}")
    return ok


def verify(sat_z_max=12, brute_max=12, log_path=DEFAULT_LOG,
           n24_sat_cap=None):
    """Cold pilot run: brute + SAT for N=1..12 (both paths), plus N=24
    via numpy brute AND SAT. All three counts must agree at N=24.

    Records every check in a deterministic JSON log (SHA-256 of the
    non-timing content is byte-stable across cold reruns).
    """
    results = {
        "task": "TASK-84 pilot: Butson N-th root vanishing-sum enumerator",
        "reference": "MUBs/BUTSON_NOTES.md",
        "reference_counts": {"Z": dict(REF_Z), "M": dict(REF_M)},
        "brute_pure": {},
        "brute_numpy": {},
        "sat_Z": {},
        "n24": {},
    }
    for n in range(1, brute_max + 1):
        zb, mb = brute_counts(n)
        results["brute_pure"][str(n)] = {"Z": zb, "M": mb}
        zn, mn = brute_counts_numpy(n)
        results["brute_numpy"][str(n)] = {"Z": zn, "M": mn}
        assert zb == zn and mb == mn, \
            f"pure vs numpy brute disagreement at N={n}: pure=({zb},{mb}) numpy=({zn},{mn})"
        if n in REF_Z:
            assert zb == REF_Z[n], f"brute Z_{n} mismatch: {zb} vs {REF_Z[n]}"
        if n in REF_M:
            assert mb == REF_M[n], f"brute M_{n} mismatch: {mb} vs {REF_M[n]}"

    for n in range(1, sat_z_max + 1):
        if n == 1:
            results["sat_Z"][str(n)] = {"Z": 0, "note": "trivial: no tail vars"}
            continue
        zs, capped, dt = sat_enumerate_Z(n)
        results["sat_Z"][str(n)] = {"Z": zs, "capped": capped}
        assert not capped
        assert zs == REF_Z[n], f"SAT Z_{n} mismatch: {zs} vs {REF_Z[n]}"

    # N = 24: numpy brute counts both Z and M; SAT independently counts Z.
    n = 24
    t0 = time.perf_counter()
    zn24, mn24 = brute_counts_numpy(n)
    results["n24"]["brute_numpy"] = {"Z": zn24, "M": mn24,
                                     "seconds": round(time.perf_counter() - t0, 3)}

    t0 = time.perf_counter()
    zs24, capped24, dt24 = sat_enumerate_Z(n, cap=n24_sat_cap)
    results["n24"]["sat_Z"] = {"Z": zs24, "capped": capped24,
                               "seconds": round(dt24, 3),
                               "cap": n24_sat_cap}
    if not capped24:
        assert zs24 == zn24, f"N=24 SAT/numpy Z mismatch: {zs24} vs {zn24}"
        results["n24"]["cross_check"] = "SAT_Z equals brute_numpy_Z"
    else:
        results["n24"]["cross_check"] = "SAT capped; no cross-check"

    # Deterministic content SHA-256: exclude wall-clock timings.
    payload = _strip_timings(results)
    payload_bytes = json.dumps(payload, sort_keys=True,
                               separators=(",", ":")).encode()
    results["deterministic_sha256"] = hashlib.sha256(payload_bytes).hexdigest()

    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)
    print(f"wrote {log_path}")
    print(f"deterministic SHA-256: {results['deterministic_sha256']}")
    return results


def _strip_timings(results):
    """Return a copy of results with all *_seconds fields removed."""
    def clean(obj):
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items()
                    if not (isinstance(k, str) and (k.endswith("seconds") or k == "elapsed"))}
        if isinstance(obj, list):
            return [clean(v) for v in obj]
        return obj
    return clean(results)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--brute", type=int, metavar="N")
    ap.add_argument("--sat-z", type=int, metavar="N")
    ap.add_argument("--log-path", default=DEFAULT_LOG)
    ap.add_argument("--sat-z-max", type=int, default=12)
    ap.add_argument("--brute-max", type=int, default=12)
    ap.add_argument("--n24-sat-cap", type=int, default=None,
                    help="cap N=24 SAT enumeration to this many solutions")
    args = ap.parse_args()

    if args.self_test:
        ok = self_test()
        sys.exit(0 if ok else 1)

    if args.brute is not None:
        zc, mc = brute_counts(args.brute)
        print(f"N={args.brute}: |Z|={zc}, |M|={mc}")
        return

    if args.sat_z is not None:
        zs, capped, dt = sat_enumerate_Z(args.sat_z)
        print(f"N={args.sat_z}: SAT |Z|={zs} (capped={capped}, {dt:.2f}s)")
        return

    if args.verify:
        verify(sat_z_max=args.sat_z_max, brute_max=args.brute_max,
               log_path=args.log_path,
               n24_sat_cap=args.n24_sat_cap)
        return

    ap.print_help()


if __name__ == "__main__":
    main()

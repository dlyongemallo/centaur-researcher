#!/usr/bin/env python3
"""Independent exact N=24 Butson quadruple census.

This file deliberately does not import either N=24 search or referee
implementation.  Arithmetic is in the power basis of Q(zeta_24), with
Phi_24 = x^8 - x^4 + 1, and graph cliques use a local bit-set search.
"""
import hashlib
import json
import os
import sys
import time
import numpy as np

N = 24
D = 6
OUT = os.path.join(os.path.dirname(__file__), "../logs/butson_n24_independent.json")

def reduction_table():
    # Multiply by x and reduce x^8 = x^4 - 1.
    table = []
    cur = [1] + [0] * 7
    for _ in range(N):
        table.append(tuple(cur))
        nxt = [0] * 8
        for j in range(7):
            nxt[j + 1] = cur[j]
        top = cur[7]
        nxt[4] += top
        nxt[0] -= top
        cur = nxt
    return np.asarray(table, dtype=np.int16)

def packed(a):
    # 24^6 < 2^28, so this is injective in signed int64 arithmetic.
    return np.asarray(a, dtype=np.int64) @ (N ** np.arange(D, dtype=np.int64))

def ray_classes(R):
    # Enumerate the dephased 24^5 tuples in slabs, exact in Z[zeta_24].
    tail = np.stack(np.meshgrid(*([np.arange(N)] * 4), indexing="ij"), -1)
    tail = tail.reshape(-1, 4).astype(np.int16)
    z = []
    m = []
    for first in range(N):
        a = np.concatenate((np.zeros((len(tail), 1), np.int16),
                            np.full((len(tail), 1), first, np.int16), tail), 1)
        coeff = R[a].sum(1)
        hit_z = np.all(coeff == 0, axis=1)
        # |sum z^a|^2 = sum_{i,j} z^(a_i-a_j).
        delta = (a[:, :, None] - a[:, None, :]) % N
        norm = R[delta].sum((1, 2))
        hit_m = (norm[:, 0] == D) & np.all(norm[:, 1:] == 0, axis=1)
        z.extend(tuple(int(x) for x in row) for row in a[hit_z])
        m.extend(tuple(int(x) for x in row) for row in a[hit_m])
    return z, m

def graph(items, zkeys):
    arr = np.asarray(items, dtype=np.int16)
    keys = np.sort(np.asarray(list(zkeys), dtype=np.int64))
    adj = [0] * len(arr)
    for lo in range(0, len(arr), 256):
        hi = min(lo + 256, len(arr))
        d = (arr[lo:hi, None, :] - arr[None, :, :]) % N
        k = packed(d)
        pos = np.searchsorted(keys, k).clip(0, len(keys) - 1)
        ok = keys[pos] == k
        for local in range(hi - lo):
            i = lo + local
            js = np.flatnonzero(ok[local, i + 1:]) + i + 1
            bits = 0
            for j in js:
                bits |= 1 << int(j)
            adj[i] = bits
    return adj

def cliques(adj, size):
    out = []
    n = len(adj)
    allbits = (1 << n) - 1
    def rec(chosen, cand):
        if len(chosen) == size:
            out.append(tuple(chosen)); return
        need = size - len(chosen)
        if cand.bit_count() < need:
            return
        while cand:
            bit = cand & -cand
            cand ^= bit
            v = bit.bit_length() - 1
            rec(chosen + [v], cand & adj[v])
    rec([], allbits)
    return out

def main():
    started = time.time()
    R = reduction_table()
    Z, M = ray_classes(R)
    zset = set(Z); mset = set(M)
    assert len(Z) == 8350 and len(M) == 19320
    zkeys = {int(packed([r])[0]) for r in Z}
    row_adj = graph(Z, zkeys)
    h_rows = cliques(row_adj, 5)
    assert len(h_rows) == 7584
    zero = (0,) * D
    census = {}
    total = 0
    max_u = max_b = 0
    rejected_pairs = 0
    for count, hs in enumerate(h_rows):
        rows = [zero] + [Z[i] for i in hs]
        cols = list(zip(*rows))
        U = []
        for r in M:
            if all(tuple((r[j] - cols[j][k]) % N for j in range(D)) in mset
                   for k in range(1, D)):
                U.append(r)
        u_adj = graph(U, zkeys)
        bases = cliques(u_adj, D)
        max_u = max(max_u, len(U)); max_b = max(max_b, len(bases))
        total += len(bases)
        census[f"{len(U)}/{len(bases)}"] = census.get(f"{len(U)}/{len(bases)}", 0) + 1
        # A candidate pair must have every cross-column difference in M.
        for ia in range(len(bases)):
            for ib in range(ia + 1, len(bases)):
                good = True
                for ca in bases[ia]:
                    for cb in bases[ib]:
                        if tuple((U[ca][j] - U[cb][j]) % N for j in range(D)) not in mset:
                            good = False; break
                    if not good: break
                if good:
                    raise AssertionError("independent pipeline found a quadruple")
                rejected_pairs += 1
    result = {
        "N": N, "zero_rays": len(Z), "mu_rays": len(M),
        "hadamards": len(h_rows), "total_completions": total,
        "max_u": max_u, "max_completions": max_b,
        "types": dict(sorted(census.items())),
        "candidate_pairs_rejected": rejected_pairs,
        "verdict": "no_quadruple", "reason": "all_completion_pairs_excluded",
        "seconds": round(time.time() - started, 3),
    }
    stable = dict(result)
    stable.pop("seconds", None)
    result["deterministic_content_sha256"] = hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
    print(json.dumps(result, indent=2, sort_keys=True))

def verify():
    with open(OUT) as f: d = json.load(f)
    assert d["N"] == 24 and d["verdict"] == "no_quadruple"
    assert d["zero_rays"] == 8350 and d["mu_rays"] == 19320
    assert d["hadamards"] == 7584 and d["total_completions"] == 1920
    assert d["max_u"] == 40 and d["max_completions"] == 4
    assert d["types"] == {"0/0": 6852, "12/4": 60, "24/4": 360, "40/0": 72, "6/1": 240}
    assert d["candidate_pairs_rejected"] == 2520
    stable = dict(d); stable.pop("seconds", None); claimed = stable.pop("deterministic_content_sha256")
    assert hashlib.sha256(json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()).hexdigest() == claimed
    print("VERIFIED independent N=24 certificate")

if __name__ == "__main__":
    if "--verify" in sys.argv: verify()
    else: main()

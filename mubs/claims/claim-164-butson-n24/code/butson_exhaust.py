"""Exact exhaustion of Butson MUB quadruples in dimension six.

For a fixed N, a BH(6,N) matrix has entries among the N-th roots of unity.
This program decides whether there is a quadruple

    {I, H1/sqrt(6), H2/sqrt(6), H3/sqrt(6)}

with all entries of H1, H2, H3 among those roots.  It uses only exact
arithmetic in Z[zeta_N], through mub_common.Cyclo.

Gauge and finite reduction
--------------------------
After root-of-unity row and column rephasings, H1 may be dephased.  Write
each root-valued ray as an exponent vector v in (Z/NZ)^6 and set v[0] = 0
to remove its scalar phase.  Define

    Z = {v : sum_j zeta_N**v[j] = 0},
    M = {v : |sum_j zeta_N**v[j]|^2 = 6}.

The five noninitial rows of H1 form a 5-clique in the graph on Z in which
v and w are adjacent iff w-v is in Z.  For each such H1, the rays unbiased
to both I and H1 are exactly

    U(H1) = {v in M : v-column(H1,c) is in M for every c}.

Orthonormal third bases are 6-cliques in U(H1), again using differences in
Z.  Two such cliques form a quadruple precisely when every cross-difference
is in M.  Thus the search below is exhaustive; row ordering is quotiented
by sorting the five rows, which is a common coordinate permutation and
does not affect existence.

Default run:

    python3 butson_exhaust.py

On the 2026-07-24 reference machine this checks N = 1,...,12 in about
twenty seconds.  The nontrivial orders are N=8 and N=12.
"""

import argparse
import itertools
import json
import time

import mub_common as mc


D = 6

# Regression values from the first independent exact run.  The "types" map
# records (number of rays in U(H1), number of six-cliques) -> number of H1.
EXPECTED = {
    1: {"zero_rays": 0, "mu_rays": 0},
    2: {"zero_rays": 10, "mu_rays": 0},
    3: {"zero_rays": 30, "mu_rays": 0},
    4: {"zero_rays": 100, "mu_rays": 0},
    5: {"zero_rays": 0, "mu_rays": 0},
    6: {"zero_rays": 340, "mu_rays": 0},
    7: {"zero_rays": 0, "mu_rays": 0},
    8: {
        "zero_rays": 640,
        "mu_rays": 1680,
        "hadamards": 432,
        "total_completions": 0,
        "max_u": 8,
        "max_completions": 0,
        "types": {"0/0": 252, "8/0": 180},
    },
    9: {"zero_rays": 270, "mu_rays": 0},
    10: {"zero_rays": 1090, "mu_rays": 0},
    11: {"zero_rays": 0, "mu_rays": 0},
    12: {
        "zero_rays": 1930,
        "mu_rays": 3600,
        "hadamards": 2184,
        "total_completions": 480,
        "max_u": 12,
        "max_completions": 4,
        "types": {"0/0": 1884, "6/1": 240, "12/4": 60},
    },
}


def ray_classes(n):
    """Return the exact sets Z and M of dephased exponent vectors."""
    roots = [mc.Cyclo.root(n, k) for k in range(n)]
    six = mc.Cyclo.integer(n, D)
    zero_rays = []
    mu_rays = []
    for tail in itertools.product(range(n), repeat=D - 1):
        vector = (0,) + tail
        value = mc.Cyclo(n)
        for exponent in vector:
            value = value + roots[exponent]
        if value.is_zero():
            zero_rays.append(vector)
        if mc.cyclo_equal(value * value.conj(), six):
            mu_rays.append(vector)
    return zero_rays, mu_rays


def difference(v, w, n):
    """w-v in (Z/nZ)^6."""
    return tuple((b - a) % n for a, b in zip(v, w))


def forward_graph(vertices, allowed_differences, n):
    """Forward adjacency bitsets; each edge is stored at its smaller end."""
    adjacency = [0] * len(vertices)
    for i, v in enumerate(vertices):
        bits = 0
        for j in range(i + 1, len(vertices)):
            if difference(v, vertices[j], n) in allowed_differences:
                bits |= 1 << j
        adjacency[i] = bits
    return adjacency


def cliques(adjacency, size):
    """Enumerate increasing cliques of a fixed size."""
    answer = []

    def extend(chosen, candidates, remaining):
        if remaining == 0:
            answer.append(tuple(chosen))
            return
        while candidates.bit_count() >= remaining:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            extend(
                chosen + [vertex],
                candidates & adjacency[vertex],
                remaining - 1,
            )

    extend([], (1 << len(adjacency)) - 1, size)
    return answer


def completion_rays(h_rows, mu_rays, mu_set, n):
    """The exact biunimodular root rays U(H)."""
    columns = list(zip(*h_rows))
    return [
        v
        for v in mu_rays
        if all(difference(column, v, n) in mu_set for column in columns)
    ]


def compatible_completions(first, second, rays, mu_set, n):
    """Whether two orthonormal ray cliques are mutually unbiased."""
    return all(
        difference(rays[i], rays[j], n) in mu_set
        for i in first
        for j in second
    )


def exact_hit_check(n, h_rows, rays, first, second):
    """Independent exact matrix-level verification of a reported hit."""
    h1 = mc.root_matrix(n, h_rows)
    h2_exp = [[rays[j][i] for j in first] for i in range(D)]
    h3_exp = [[rays[j][i] for j in second] for i in range(D)]
    h2 = mc.root_matrix(n, h2_exp)
    h3 = mc.root_matrix(n, h3_exp)
    assert all(mc.exact_orthonormal(h, D) for h in (h1, h2, h3))
    assert all(
        mc.exact_mu_pair(a, b, D)
        for a, b in ((h1, h2), (h1, h3), (h2, h3))
    )
    return {"H1": h_rows, "H2": h2_exp, "H3": h3_exp}


def exhaust_order(n, progress=False):
    """Return an exact verdict and census for one root order."""
    started = time.time()
    zero_rays, mu_rays = ray_classes(n)
    result = {
        "N": n,
        "zero_rays": len(zero_rays),
        "mu_rays": len(mu_rays),
        "verdict": "no_quadruple",
    }

    # If Z is empty there is no BH(6,N).  If M is empty, no two
    # root-valued bases can be mutually unbiased.
    if not zero_rays or not mu_rays:
        result["reason"] = (
            "no_BH(6,N)" if not zero_rays else "no_unbiased_root_sum"
        )
        result["seconds"] = round(time.time() - started, 3)
        check_regression(result)
        return result

    zero_set = set(zero_rays)
    mu_set = set(mu_rays)
    row_graph = forward_graph(zero_rays, zero_set, n)
    hadamard_row_sets = cliques(row_graph, D - 1)
    result["hadamards"] = len(hadamard_row_sets)

    census = {}
    total_completions = 0
    max_u = 0
    max_completions = 0
    zero_row = (0,) * D

    for h_index, row_set in enumerate(hadamard_row_sets):
        h_rows = [zero_row] + [zero_rays[i] for i in row_set]
        rays = completion_rays(h_rows, mu_rays, mu_set, n)
        max_u = max(max_u, len(rays))
        ray_graph = forward_graph(rays, zero_set, n)
        bases = cliques(ray_graph, D)
        max_completions = max(max_completions, len(bases))
        total_completions += len(bases)
        key = "%d/%d" % (len(rays), len(bases))
        census[key] = census.get(key, 0) + 1

        for a in range(len(bases)):
            for b in range(a + 1, len(bases)):
                if compatible_completions(
                    bases[a], bases[b], rays, mu_set, n
                ):
                    result["verdict"] = "quadruple"
                    result["witness"] = exact_hit_check(
                        n, h_rows, rays, bases[a], bases[b]
                    )
                    result["seconds"] = round(time.time() - started, 3)
                    return result

        if progress and (h_index + 1) % 500 == 0:
            print(
                "N=%d: checked %d/%d H1 candidates"
                % (n, h_index + 1, len(hadamard_row_sets)),
                flush=True,
            )

    result.update(
        {
            "total_completions": total_completions,
            "max_u": max_u,
            "max_completions": max_completions,
            "types": dict(sorted(census.items())),
            "reason": "all_completion_pairs_excluded",
            "seconds": round(time.time() - started, 3),
        }
    )
    check_regression(result)
    return result


def check_regression(result):
    """Fail loudly if a reference census changes."""
    expected = EXPECTED.get(result["N"])
    if expected is None:
        return
    for key, value in expected.items():
        if result.get(key) != value:
            raise AssertionError(
                "N=%d regression for %s: got %r, expected %r"
                % (result["N"], key, result.get(key), value)
            )


def parse_orders(text):
    orders = []
    for part in text.split(","):
        if "-" in part:
            lo, hi = (int(x) for x in part.split("-", 1))
            orders.extend(range(lo, hi + 1))
        else:
            orders.append(int(part))
    if any(n < 1 for n in orders):
        raise ValueError("root orders must be positive")
    return sorted(set(orders))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orders",
        default="1-12",
        help="comma-separated root orders and ranges (default: 1-12)",
    )
    parser.add_argument("--json", action="store_true", help="emit JSON lines")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()

    results = []
    for n in parse_orders(args.orders):
        result = exhaust_order(n, progress=args.progress)
        results.append(result)
        if args.json:
            print(json.dumps(result, sort_keys=True))
        else:
            print(
                "N=%d: %s; |Z|=%d, |M|=%d, H=%s, completions=%s (%.3fs)"
                % (
                    n,
                    result["verdict"],
                    result["zero_rays"],
                    result["mu_rays"],
                    result.get("hadamards", "-"),
                    result.get("total_completions", "-"),
                    result["seconds"],
                )
            )

    if not args.json:
        if all(r["verdict"] == "no_quadruple" for r in results):
            print("CERTIFIED: no requested Butson MUB quadruple exists.")
        else:
            print("EXACT COUNTEREXAMPLE FOUND.")


if __name__ == "__main__":
    main()

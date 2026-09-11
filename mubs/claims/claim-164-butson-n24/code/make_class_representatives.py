"""Generate data/bh6_24_class_representatives.txt.

Enumerates the 7,584 dephased BH(6, 24) row sets, partitions them into
monomial-equivalence classes via the moves in
`bh6_full_monomial_classes.py`, and writes one block per class: the
lexicographically minimal member as a 6 x 6 exponent matrix over Z/24
(row 0 is the all-ones row), together with the class size, the
completion census (|U(H)| and number of completion bases), the proper
divisors m of 24 at which the class has a representative in BH(6, m),
and the class of the transposed representative. Output is
deterministic: classes are sorted by their minimal member.
"""

import os
from collections import deque

from bh6_full_monomial_classes import neighbours
from butson_bh6_equivalence_classes import (
    D, enumerate_dephased_hadamards, transpose_rowset)
from butson_n24_search import ray_classes_numpy
import butson_exhaust as be

N = 24
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "../data/bh6_24_class_representatives.txt")


def main():
    zero_rays, rowsets, _ = enumerate_dephased_hadamards(N)
    _, mu_rays = ray_classes_numpy(N)
    zero_set, mu_set = set(zero_rays), set(mu_rays)
    mu_rays = sorted(mu_rays)

    index = {rs: k for k, rs in enumerate(rowsets)}
    comp = [-1] * len(rowsets)
    ncomp = 0
    for k, rs in enumerate(rowsets):
        if comp[k] >= 0:
            continue
        comp[k] = ncomp
        queue = deque([rs])
        while queue:
            cur = queue.popleft()
            for nb in neighbours(cur, N, False, False):
                j = index[nb]
                if comp[j] < 0:
                    comp[j] = ncomp
                    queue.append(nb)
        ncomp += 1

    members = [[] for _ in range(ncomp)]
    for k, rs in enumerate(rowsets):
        members[comp[k]].append(rs)
    # Deterministic order and ids: sort classes by minimal member.
    order = sorted(range(ncomp), key=lambda c: min(members[c]))
    new_id = {c: i for i, c in enumerate(order)}

    def profile(rep):
        h_rows = [(0,) * D] + list(rep)
        u = be.completion_rays(h_rows, mu_rays, mu_set, N)
        if not u:
            return 0, 0
        bases = be.cliques(be.forward_graph(sorted(u), zero_set, N), D)
        return len(u), len(bases)

    def visible(cls):
        out = []
        for m in (2, 3, 4, 6, 8, 12):
            step = N // m
            if any(all(x % step == 0 for v in rs for x in v)
                   for rs in members[cls]):
                out.append(m)
        return out

    lines = [
        "# The 25 monomial-equivalence classes of BH(6, 24).",
        "# One block per class: 'class <id>: size=<n> U=<u> bases=<b> "
        "sub_orders=<m,...> transpose_class=<id>'",
        "# followed by the lexicographically minimal dephased member as a "
        "6 x 6 matrix of",
        "# exponents e (entry = exp(2*pi*i*e/24)); row 0 and column 0 are "
        "the dephasing anchors.",
        "# Generated deterministically by code/make_class_representatives.py.",
        "",
    ]
    for c in order:
        rep = min(members[c])
        usz, nb = profile(rep)
        vis = visible(c) or ["-"]
        t_cls = new_id[comp[index[transpose_rowset(rep)]]]
        lines.append(
            f"class {new_id[c]}: size={len(members[c])} U={usz} bases={nb} "
            f"sub_orders={','.join(str(m) for m in vis)} "
            f"transpose_class={t_cls}")
        for row in ((0,) * D,) + rep:
            lines.append(" ".join(f"{x:2d}" for x in row))
        lines.append("")
    with open(OUT, "w") as f:
        f.write("\n".join(lines))
    print(f"wrote {OUT}: {ncomp} classes")


if __name__ == "__main__":
    main()

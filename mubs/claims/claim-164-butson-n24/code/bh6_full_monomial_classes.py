"""Full monomial-equivalence class counts for BH(6, N) dephased row sets.

Literature-comparison driver for TASK-111 (host-side, network-era follow-up).
Extends `butson_bh6_equivalence_classes.py`'s residual group (S5 column
permutations, optional transpose/conjugation) to the FULL residual group of
monomial equivalence acting on dephased representatives, by adding the two
re-dephasing moves:

  - column-anchor move: swap coordinate 0 with coordinate 1, then re-dephase
    every row by subtracting its new 0th entry (row rephasing);
  - row-anchor move: make row r the all-ones row, i.e., replace the set S by
    {-r} + {v - r : v in S \\ {r}} (column rephasing by zeta^{-r_j}).

Together with S5 on coordinates 1..5 and the free row-sort, these generate
the full orbit of dephased forms under monomial equivalence with N-th-root
monomial matrices (Lampio-Ostergard-Szollosi convention).  Reported per N:
class counts under M (monomial), M+T (transpose), M+K (conjugation), M+T+K.
"""

import sys
from collections import deque

from butson_bh6_equivalence_classes import (
    D,
    enumerate_dephased_hadamards,
    transpose_rowset,
    conjugate_rowset,
)


def col_swap_adj(rows, i):
    """Transpose coordinates i and i+1 (1 <= i <= 4). Stays dephased."""
    new_rows = []
    for v in rows:
        w = list(v)
        w[i], w[i + 1] = w[i + 1], w[i]
        new_rows.append(tuple(w))
    return tuple(sorted(new_rows))


def col_anchor_move(rows, n):
    """Swap coordinates 0 and 1, then re-dephase each row."""
    new_rows = []
    for v in rows:
        w = (v[1], 0) + v[2:]
        s = w[0]
        new_rows.append(tuple((x - s) % n for x in w))
    return tuple(sorted(new_rows))


def row_anchor_move(rows, r, n):
    """Make row r the all-ones row: S -> {-r} + {v - r : v != r}."""
    new_rows = [tuple((-x) % n for x in r)]
    for v in rows:
        if v == r:
            continue
        new_rows.append(tuple((a - b) % n for a, b in zip(v, r)))
    return tuple(sorted(new_rows))


def neighbours(rows, n, include_T, include_K):
    for i in range(1, D - 1):
        yield col_swap_adj(rows, i)
    yield col_anchor_move(rows, n)
    for r in rows:
        yield row_anchor_move(rows, r, n)
    if include_T:
        yield transpose_rowset(rows)
    if include_K:
        yield conjugate_rowset(rows, n)


def class_count(rowsets, n, include_T, include_K):
    index = {rs: k for k, rs in enumerate(rowsets)}
    seen = [False] * len(rowsets)
    classes = 0
    for k, rs in enumerate(rowsets):
        if seen[k]:
            continue
        classes += 1
        queue = deque([rs])
        seen[k] = True
        while queue:
            cur = queue.popleft()
            for nb in neighbours(cur, n, include_T, include_K):
                j = index.get(nb)
                if j is None:
                    raise AssertionError(
                        f"orbit left the enumerated set at N={n}: {nb}")
                if not seen[j]:
                    seen[j] = True
                    queue.append(nb)
    return classes


def main():
    orders = [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1
                               else "3,4,6,8,12,24".split(","))]
    print("N | dephased | M | M+T | M+K | M+T+K")
    for n in orders:
        _, rowsets, _ = enumerate_dephased_hadamards(n)
        row = [str(n), str(len(rowsets))]
        for include_T, include_K in ((False, False), (True, False),
                                     (False, True), (True, True)):
            row.append(str(class_count(rowsets, n, include_T, include_K)))
        print(" | ".join(row))


if __name__ == "__main__":
    main()

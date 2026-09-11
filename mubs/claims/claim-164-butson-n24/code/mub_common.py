"""Core toolkit for the d = 6 MUB project (stdlib only).

Conventions. A basis of C^d is represented by its "unimodular representative":
a d x d complex matrix B with all entries of modulus 1 such that B B* = d I,
i.e., B / sqrt(d) is unitary (a rescaled complex Hadamard matrix, CHM). The
computational basis is represented by sqrt(d) X = B with B = "the identity
basis" handled separately, since its representative is not unimodular; use
`Basis` below, which stores the unitary directly.

Two unitaries X, Y are mutually unbiased iff sqrt(d) X* Y is a CHM with
unimodular entries. For triplets we follow Matolcsi-Matszangosz-Varga-Weiner
(arXiv:2503.14752, "MMVW"): the transition matrices of a triplet (X1, X2, X3)
are H1 = sqrt(d) X1* X2, H2 = sqrt(d) X2* X3, H3 = sqrt(d) X3* X1.

The exact layer works in Z[zeta_N] (N-th roots of unity) with integer
coefficients, tracking powers of sqrt(6) separately, so that MU properties and
the MMVW identity can be verified exactly (certified) for matrices whose
entries are roots of unity over sqrt(6).
"""

import cmath
import itertools
import math

TOL = 1e-9


# ---------------------------------------------------------------------------
# Float complex matrices (lists of lists).
# ---------------------------------------------------------------------------

def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    Bt = list(zip(*B))
    return [[sum(A[i][k] * Bt[j][k] for k in range(m)) for j in range(p)]
            for i in range(n)]


def dagger(A):
    return [[A[j][i].conjugate() for j in range(len(A))]
            for i in range(len(A[0]))]


def eye(n):
    return [[1.0 + 0j if i == j else 0j for j in range(n)] for i in range(n)]


def mat_scale(A, s):
    return [[s * x for x in row] for row in A]


def mat_max_dev(A, B):
    return max(abs(A[i][j] - B[i][j])
               for i in range(len(A)) for j in range(len(A[0])))


def is_chm(H, tol=TOL):
    """Entries unimodular and H H* = d I (rescaled complex Hadamard)."""
    d = len(H)
    if any(abs(abs(x) - 1.0) > tol for row in H for x in row):
        return False
    return mat_max_dev(mat_mul(H, dagger(H)), mat_scale(eye(d), d)) < d * tol


def dephase(H):
    """Multiply rows and columns by phases so row 0 and column 0 become 1."""
    d = len(H)
    out = [[H[i][j] / H[i][0] / H[0][j] * H[0][0] for j in range(d)]
           for i in range(d)]
    return out


# ---------------------------------------------------------------------------
# MU checks and the MMVW loss functional.
# ---------------------------------------------------------------------------

def mu_pair_deficit(X, Y):
    """Sum of squared deviations of |(X* Y)_kl| from 1/sqrt(d)."""
    d = len(X)
    M = mat_mul(dagger(X), Y)
    t = 1.0 / math.sqrt(d)
    return sum((abs(M[k][l]) - t) ** 2 for k in range(d) for l in range(d))


def orthonormality_deficit(X):
    d = len(X)
    M = mat_mul(dagger(X), X)
    return sum(abs(M[k][l] - (1.0 if k == l else 0.0)) ** 2
               for k in range(d) for l in range(d))


def mub_loss(bases):
    """The MMVW loss: zero iff the unitaries form a set of MUBs.

    bases is a list of d x d matrices interpreted as unitaries (columns are
    the basis vectors).
    """
    loss = sum(orthonormality_deficit(X) for X in bases)
    for X, Y in itertools.combinations(bases, 2):
        loss += mu_pair_deficit(X, Y)
    return loss


def constellation_deficit(blocks, d=None):
    """Brierley-Weigert style deficit for an MU constellation.

    blocks is a list of lists of vectors (tuples of complex numbers). Vectors
    within a block should be orthonormal; vectors in different blocks should
    have |<u,v>|^2 = 1/d. Returns the sum of squared constraint violations.
    """
    if d is None:
        d = len(blocks[0][0])
    t = 1.0 / d
    total = 0.0
    for bi, block in enumerate(blocks):
        for vi, v in enumerate(block):
            nrm = sum(abs(x) ** 2 for x in v)
            total += (nrm - 1.0) ** 2
            for w in block[vi + 1:]:
                ip = sum(a.conjugate() * b for a, b in zip(v, w))
                total += abs(ip) ** 2
            for block2 in blocks[bi + 1:]:
                for w in block2:
                    ip = sum(a.conjugate() * b for a, b in zip(v, w))
                    total += (abs(ip) ** 2 - t) ** 2
    return total


# ---------------------------------------------------------------------------
# Karlsson H2-reducibility test and the Haagerup invariant.
# ---------------------------------------------------------------------------

def is_h2_reducible(H, tol=1e-7):
    """Karlsson (arXiv:1003.4133): a CHM of order 6 is H2-reducible iff it
    contains a 2x2 Hadamard submatrix, i.e., rows a < b and columns i < j
    with z_ai z_bj + z_aj z_bi = 0."""
    d = len(H)
    for a in range(d):
        for b in range(a + 1, d):
            for i in range(d):
                for j in range(i + 1, d):
                    if abs(H[a][i] * H[b][j] + H[a][j] * H[b][i]) < tol:
                        return True
    return False


def haagerup_invariant(H, ndigits=6):
    """The multiset of phases of h_ij h_kl conj(h_il) conj(h_kj), rounded.

    Invariant under the CHM equivalence relation (rephasing and permutations
    of rows and columns, and global phase); a fingerprint for orbit testing.
    Necessary, not sufficient: equal fingerprints do not prove equivalence.
    """
    d = len(H)
    vals = []
    for i in range(d):
        for k in range(d):
            for j in range(d):
                for l in range(d):
                    z = (H[i][j] * H[k][l] *
                         H[i][l].conjugate() * H[k][j].conjugate())
                    vals.append((round(z.real, ndigits), round(z.imag, ndigits)))
    return tuple(sorted(vals))


# ---------------------------------------------------------------------------
# Exact layer: Z[zeta_N] with integer coefficients.
#
# An element is a tuple of N integers (c_0, ..., c_{N-1}) representing
# sum c_k zeta_N^k, i.e., an integer polynomial modulo x^N - 1. Zero testing
# reduces modulo the N-th cyclotomic polynomial (monic with integer
# coefficients, so the reduction stays over the integers).
# ---------------------------------------------------------------------------

def _poly_divmod_int(num, den):
    """Divide integer polynomials (lists, low degree first); den monic."""
    num = list(num)
    dd = len(den) - 1
    while len(den) > 1 and den[-1] == 0:
        den = den[:-1]
        dd -= 1
    q = [0] * max(1, len(num) - dd)
    while len(num) - 1 >= dd and any(num):
        while num and num[-1] == 0:
            num.pop()
        if len(num) - 1 < dd:
            break
        shift = len(num) - 1 - dd
        c = num[-1]
        q[shift] += c
        for i, dc in enumerate(den):
            num[shift + i] -= c * dc
    while num and num[-1] == 0:
        num.pop()
    return q, num


_CYCLO_CACHE = {}


def cyclotomic_poly(n):
    """Integer coefficient list (low degree first) of Phi_n(x)."""
    if n in _CYCLO_CACHE:
        return _CYCLO_CACHE[n]
    # x^n - 1 divided by the product of Phi_d for proper divisors d of n.
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            q, r = _poly_divmod_int(num, cyclotomic_poly(d))
            assert not r, "cyclotomic division must be exact"
            num = q
    _CYCLO_CACHE[n] = num
    return num


class Cyclo:
    """Element of Z[zeta_N], stored modulo x^N - 1 with integer coeffs."""

    __slots__ = ("N", "c")

    def __init__(self, N, coeffs=None):
        self.N = N
        self.c = [0] * N if coeffs is None else list(coeffs)

    @classmethod
    def root(cls, N, k):
        z = cls(N)
        z.c[k % N] = 1
        return z

    @classmethod
    def integer(cls, N, n):
        z = cls(N)
        z.c[0] = n
        return z

    def __add__(self, other):
        return Cyclo(self.N, [a + b for a, b in zip(self.c, other.c)])

    def __sub__(self, other):
        return Cyclo(self.N, [a - b for a, b in zip(self.c, other.c)])

    def __mul__(self, other):
        if isinstance(other, int):
            return Cyclo(self.N, [a * other for a in self.c])
        N = self.N
        out = [0] * N
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(other.c):
                    if b:
                        out[(i + j) % N] += a * b
        return Cyclo(N, out)

    def conj(self):
        N = self.N
        out = [0] * N
        for i, a in enumerate(self.c):
            out[(-i) % N] += a
        return Cyclo(N, out)

    def is_zero(self):
        if not any(self.c):
            return True
        _, r = _poly_divmod_int(self.c, cyclotomic_poly(self.N))
        return not r

    def to_complex(self):
        return sum(a * cmath.exp(2j * math.pi * k / self.N)
                   for k, a in enumerate(self.c) if a)

    def __repr__(self):
        return "Cyclo(%d, %r)" % (self.N, self.c)


def cyclo_equal(a, b):
    return (a - b).is_zero()


# ---------------------------------------------------------------------------
# Exact matrices over Z[zeta_N]: entries are Cyclo values. A "root matrix"
# has entries that are single roots of unity.
# ---------------------------------------------------------------------------

def root_matrix(N, exponents):
    """Matrix of roots of unity from a matrix of integer exponents mod N."""
    return [[Cyclo.root(N, e) for e in row] for row in exponents]


def cyclo_mat_mul_dagger_left(A, B):
    """Compute A* B exactly for Cyclo matrices (A* is conjugate transpose)."""
    n = len(A)
    m = len(A[0])
    p = len(B[0])
    out = []
    for i in range(m):
        row = []
        for j in range(p):
            acc = Cyclo(A[0][0].N)
            for k in range(n):
                acc = acc + A[k][i].conj() * B[k][j]
            row.append(acc)
        out.append(row)
    return out


def exact_orthonormal(B, d):
    """Check B B* = d I exactly, for a Cyclo matrix with unimodular entries."""
    n = len(B)
    M = cyclo_mat_mul_dagger_left(B, B)  # B* B; for square unimodular B with
    # B* B = d I we also get B B* = d I, and B* B is what column-orthonormality
    # needs (columns are the basis vectors).
    for i in range(n):
        for j in range(n):
            target = Cyclo.integer(B[0][0].N, d if i == j else 0)
            if not cyclo_equal(M[i][j], target):
                return False
    return True


def exact_mu_pair(B1, B2, d):
    """Check |(B1* B2)_kl|^2 = d exactly (entries of B1, B2 unimodular)."""
    M = cyclo_mat_mul_dagger_left(B1, B2)
    for row in M:
        for m in row:
            if not cyclo_equal(m * m.conj(), Cyclo.integer(m.N, d)):
                return False
    return True


def cyclo_matrix_to_complex(B):
    return [[x.to_complex() for x in row] for row in B]

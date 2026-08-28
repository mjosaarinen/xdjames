"""Matrices over F_q: sampling, inversion, rank.

Matrices are held as lists of rows of small ints -- readable, and the natural
shape for the scheme's bookkeeping.  Only Gauss-Jordan runs on the packed
representation, because S and T reach 578 x 578 and an elimination written
over Python lists would be O(n^3) interpreter steps.

Packed elimination keeps lanes *unreduced* while sweeping a column: over odd
F_p a row accumulates at most n * (p-1)^2 before anything is read back, which
sits far below the 40-bit lanes used here, so no reduction is needed in the
inner loop.
"""

from .ff import F2Vec, F2kVec, FpVec
from .symmetric import sample_fq


def _packed(fq):
    if fq.q == 2:
        return F2Vec()
    if fq.p == 2:
        return F2kVec(fq)
    return FpVec(fq.q, 40)


def random_matrix(fq, nrows, ncols, xof):
    flat = sample_fq(xof, nrows * ncols, fq.q)
    return [flat[i * ncols:(i + 1) * ncols] for i in range(nrows)]


def identity(fq, n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def _echelon(fq, M, ncols, augment=None):
    """Reduced row echelon form of M (optionally carrying an augmentation).

    Returns (rank, rows, aug_rows) with rows/aug_rows as packed vectors.
    """
    V = _packed(fq)
    q = fq.q
    rows = [V.from_list(r) for r in M]
    aug = [V.from_list(r) for r in augment] if augment is not None else None
    nrows = len(rows)
    rank = 0
    for c in range(ncols):
        piv = None
        for i in range(rank, nrows):
            if V.coef(V.reduce(rows[i]), c) % q:
                piv = i
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        if aug is not None:
            aug[rank], aug[piv] = aug[piv], aug[rank]
        rows[rank] = V.reduce(rows[rank])
        pv = V.coef(rows[rank], c) % q
        if pv != 1:
            s = fq.inv(pv)
            rows[rank] = V.reduce(V.scal(rows[rank], s))
            if aug is not None:
                aug[rank] = V.reduce(V.scal(aug[rank], s))
        if aug is not None:
            aug[rank] = V.reduce(aug[rank])
        for i in range(nrows):
            if i == rank:
                continue
            v = V.coef(rows[i], c) % q
            if v:
                f = fq.neg(v)
                rows[i] = V.add(rows[i], V.scal(rows[rank], f))
                if aug is not None:
                    aug[i] = V.add(aug[i], V.scal(aug[rank], f))
        rank += 1
        if rank == nrows:
            break
    rows = [V.reduce(r) for r in rows]
    if aug is not None:
        aug = [V.reduce(r) for r in aug]
    return rank, rows, aug, V


def rank(fq, M, ncols):
    return _echelon(fq, M, ncols)[0]


def invert(fq, M):
    """Inverse of a square matrix, or None if singular."""
    n = len(M)
    rk, rows, aug, V = _echelon(fq, M, n, augment=identity(fq, n))
    if rk != n:
        return None
    return [V.to_list(r, n) for r in aug]


def random_invertible(fq, n, xof):
    """A uniform invertible n x n matrix, together with its inverse.

    Rejection sampling: a uniform matrix over F_q is invertible with
    probability prod_{i=1..n} (1 - q^-i), about 0.29 for q = 2 and rising
    quickly with q, so the loop turns over a handful of times at most.
    """
    while True:
        M = random_matrix(fq, n, n, xof)
        Minv = invert(fq, M)
        if Minv is not None:
            return M, Minv


def random_full_rank(fq, nrows, ncols, xof):
    """A uniform nrows x ncols matrix of full row rank (nrows <= ncols)."""
    while True:
        M = random_matrix(fq, nrows, ncols, xof)
        if rank(fq, M, ncols) == nrows:
            return M


def mat_vec(fq, M, v):
    """M . v  (v a column of length ncols) -> list of length nrows."""
    out = []
    for row in M:
        acc = 0
        for c, x in zip(row, v):
            if c and x:
                acc = fq.add(acc, fq.mul(c, x))
        out.append(acc)
    return out


def vec_mat(fq, v, M):
    """v . M  (v a row of length nrows) -> list of length ncols."""
    ncols = len(M[0])
    out = [0] * ncols
    for x, row in zip(v, M):
        if x:
            for j, c in enumerate(row):
                if c:
                    out[j] = fq.add(out[j], fq.mul(x, c))
    return out

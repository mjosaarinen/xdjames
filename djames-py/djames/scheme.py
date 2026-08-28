"""James and D-James: key generation, signing, verification.

Notation follows the paper.  Row-vector convention throughout: the signature
is a in F_q^n, the hidden HFE variable is x = a S, and a = x S^-1.

Central map (D-James; James is the same without the Dragon block):

    H(X, y) = sum_{(i,j) in monomials} lambda_ij X^(q^i + q^j)      HFE core
            + sum_{k <= d}   L_k(y)  X^(q^k)                        Dragon
            + sum_{k <  d}   sum_j M[k][j] z_j X^(q^k)              IP, bilinear
            + sum_{i <= j}   G[i][j] z_i z_j                        IP, quadratic

with z = M_Z x in F_q^r and L_k(y) = sum_j Lambda[k][j] y_j in K.  Every term
is quadratic in (a, b) jointly, so the public system is homogeneous -- which
is why the all-zero signature is a trivial solution and must be rejected.

Fixing the message hash y turns each L_k(y) into a constant of K, so the
central map collapses to a univariate polynomial of degree D = q^d + 1 and
inverting it costs no more than plain HFE.  That is the whole point of the
Dragon construction: the equation count m is decoupled from the hash length
n_y, so the signature can be shorter than 2*lambda bits.
"""

import itertools

from .ff import ExtField, make_vec, prep, lincomb_prep
from .linalg import random_invertible, random_full_rank, mat_vec, vec_mat
from .poly import roots, norm
from .symmetric import XOF, DOM_KEY, DOM_PRF, DOM_PAD, DOM_EDF, hash_message, sample_fq
from .codec import DigitPacker, bytes_to_digits, encode_vector, decode_vector

# The D-James signer walks salt = 0, 1, 2, ... and the verifier re-walks it,
# because the salt is not transmitted (paper, footnote 5).
#
# Treating salts as independent, the paper's expected trial counts (1.582 in
# characteristic 2, 3.07 otherwise, its footnote 3) give per-salt failure
# probabilities of 1 - 1/1.582 = 0.368 and 1 - 1/3.07 = 0.674, so 256 salts
# fail with probability about 2^-369 and 2^-146 respectively.  The odd figure
# is the one that binds.
#
# The cap is security-relevant, not merely an engineering limit: a forger may
# aim at any of MAX_SALT hash values, which concedes log2(MAX_SALT) = 8 bits
# against generic attacks.  A real specification has to pin it deliberately.
MAX_SALT = 256

_REDUCE_EVERY = 256

# A secret key is exactly its seed, so the seed carries the full security of
# the key and its length is part of the interface rather than a suggestion.
SEED_BYTES = 32


class SecretKey:
    __slots__ = ("params", "seed", "K", "lam", "Lam", "MZ", "Mbil", "G",
                 "S", "Sinv", "T", "Tinv", "prf")

    def to_bytes(self):
        """A secret key is exactly its seed: everything else is derived."""
        return self.seed


class PublicKey:
    """m homogeneous quadratic equations over F_q in (a, b).

    Coefficients are stored transposed: one packed F_q^m vector per monomial,
    holding that monomial's coefficient in all m equations at once.  Evaluation
    is then a handful of packed additions rather than m separate polynomial
    evaluations.
    """

    __slots__ = ("params", "W", "aa", "ab", "_off")

    def __init__(self, params, W, aa, ab):
        self.params, self.W, self.aa, self.ab = params, W, aa, ab
        self._off = _offsets(params.n)

    def idx_aa(self, i, j):
        if i > j:
            i, j = j, i
        return self._off[i] + j - i

    def to_bytes(self):
        P = self.params
        pk = DigitPacker(P.q)
        for v in self.aa:
            pk.push_packed(self.W, v, P.m)
        for v in self.ab:
            pk.push_packed(self.W, v, P.m)
        return pk.bytes()

    @classmethod
    def from_bytes(cls, P, data):
        W = _pk_space(P)
        total = P.pk_coeffs * P.m
        digits = bytes_to_digits(data, total, P.q)
        vecs = [W.from_list(digits[i * P.m:(i + 1) * P.m])
                for i in range(P.pk_coeffs)]
        return cls(P, W, vecs[:P.n_aa], vecs[P.n_aa:])


def _offsets(n):
    off, acc = [], 0
    for i in range(n):
        off.append(acc)
        acc += n - i
    return off


def _pk_space(P):
    from .ff import Fq
    # Lanes wide enough that evaluation can accumulate _REDUCE_EVERY terms
    # before a reduction is needed.
    return make_vec(Fq(P.q), 2 * _REDUCE_EVERY)


# --------------------------------------------------------------- key  gen


def keygen(P, seed):
    """Derive a key pair from a 32-byte seed.  Fully deterministic."""
    if len(seed) != SEED_BYTES:
        raise ValueError("seed must be exactly %d bytes, got %d"
                         % (SEED_BYTES, len(seed)))
    xof = XOF(DOM_KEY, seed, P.tag())
    K = ExtField(P.q, P.n, fpoly=P.fpoly, verify=False)
    fq = K.fq
    n, m, r, d, ny = P.n, P.m, P.r, P.d, P.ny

    sk = SecretKey()
    sk.params, sk.seed, sk.K = P, seed, K
    # Sampling order is part of the spec: changing it changes every key.
    sk.lam = [K.random_nonzero(xof) for _ in P.monomials]
    sk.Mbil = [[K.random(xof) for _ in range(r)] for _ in range(d)]
    sk.G = [[K.random(xof) if i <= j else None for j in range(r)]
            for i in range(r)]
    sk.Lam = ([[K.random(xof) for _ in range(ny)] for _ in range(d + 1)]
              if ny is not None else None)
    sk.MZ = random_full_rank(fq, r, n, xof)
    sk.S, sk.Sinv = random_invertible(fq, n, xof)
    sk.T, sk.Tinv = random_invertible(fq, n, xof)
    sk.prf = XOF(DOM_PRF, seed, P.tag()).read(32)

    return _public_key(P, sk), sk


def _linear_forms(P, sk):
    """The central map's inputs written as linear forms in (a, b) over K.

    X = sum_i a_i sigma_i with sigma_i = sum_c S[i][c] alpha^c, so X^(q^k) has
    coefficients sigma_i^(q^k) -- Frobenius is F_q-linear and the a_i lie in
    F_q.  Likewise z_j = sum_i a_i * (M_Z S^T)[j][i].
    """
    K, fq, n, r, d = sk.K, sk.K.fq, P.n, P.r, P.d
    sigma = [K.from_coords(sk.S[i]) for i in range(n)]
    A = [[K.frob(s, k) for s in sigma] for k in range(d + 1)]
    zf = []
    for j in range(r):
        row = sk.MZ[j]
        col = []
        for i in range(n):
            acc = 0
            Si = sk.S[i]
            for c in range(n):
                if row[c] and Si[c]:
                    acc = fq.add(acc, fq.mul(row[c], Si[c]))
            col.append(acc)
        zf.append(col)
    return A, zf


def _terms(P, sk, A, zf):
    """The central map as a sum of rank-one outer products over K.

    Every piece of H is a product of two F_q-linear forms, so the quadratic
    form it induces is an outer product u (x) v.  The IP families are folded
    over their inner index first (sum_j M[k][j] Z_j is itself a linear form),
    which cuts the term count from d*r + r(r+1)/2 down to d + r.
    """
    K, n, r, d = sk.K, P.n, P.r, P.d
    aa, ab = [], []

    for (mi, mj), lm in zip(P.monomials, sk.lam):
        aa.append(([K.mul(lm, x) for x in A[mi]], A[mj]))

    Z = [[K.scal(K.ONE, zf[j][i]) for i in range(n)] for j in range(r)]

    for k in range(d):
        v = [K.ZERO] * n
        for j in range(r):
            c = sk.Mbil[k][j]
            if not K.is_zero(c):
                for i in range(n):
                    if zf[j][i]:
                        v[i] = K.add(v[i], K.scal(c, zf[j][i]))
        aa.append((A[k], v))

    for i0 in range(r):
        v = [K.ZERO] * n
        for j0 in range(i0, r):
            c = sk.G[i0][j0]
            if not K.is_zero(c):
                for i in range(n):
                    if zf[j0][i]:
                        v[i] = K.add(v[i], K.scal(c, zf[j0][i]))
        aa.append((Z[i0], v))

    if P.ny is not None:
        for k in range(d + 1):
            ab.append((A[k], sk.Lam[k]))
    return aa, ab


def _public_key(P, sk):
    K, fq, V, n, m, ny = sk.K, sk.K.fq, sk.K.V, P.n, P.m, P.ny
    A, zf = _linear_forms(P, sk)
    aa_terms, ab_terms = _terms(P, sk, A, zf)
    off = _offsets(n)

    # Accumulate the K-valued coefficient of each monomial, over all terms.
    # Applying the output projection once at the end -- rather than per term --
    # roughly halves the work.
    Caa = [K.ZERO] * (n * (n + 1) // 2)
    for (u, v) in aa_terms:
        pv = [prep(V, x, n) for x in v]
        for i in range(n):
            wi = u[i]
            if K.is_zero(wi):
                continue
            tbl, cur = [wi], wi
            for _ in range(n - 1):
                cur = K.mul_t(cur)
                tbl.append(cur)
            for j in range(n):
                p = lincomb_prep(V, pv[j], tbl)
                idx = off[i] + j - i if i <= j else off[j] + i - j
                Caa[idx] = V.add(Caa[idx], p)

    Cab = [K.ZERO] * (n * ny) if ny else []
    for (u, v) in ab_terms:
        pv = [prep(V, x, n) for x in v]
        for i in range(n):
            wi = u[i]
            if K.is_zero(wi):
                continue
            tbl, cur = [wi], wi
            for _ in range(n - 1):
                cur = K.mul_t(cur)
                tbl.append(cur)
            base = i * ny
            for j in range(ny):
                Cab[base + j] = V.add(Cab[base + j], lincomb_prep(V, pv[j], tbl))

    # Project n central equations down to the m published ones: the minus
    # modifier is exactly "keep the first m rows of T".
    W = _pk_space(P)
    PI = [W.from_list([sk.T[k][c] for k in range(m)]) for c in range(n)]
    aa = [lincomb_prep(W, prep(V, V.reduce(c), n), PI) for c in Caa]
    ab = [lincomb_prep(W, prep(V, V.reduce(c), n), PI) for c in Cab]
    return PublicKey(P, W, aa, ab)


# ---------------------------------------------------------------- signing


def _central_poly(P, sk, consts, z):
    """The univariate polynomial to root-find, for one hash and one IP guess."""
    K, q, r, d = sk.K, P.q, P.r, P.d
    fq = K.fq
    F = [K.ZERO] * (P.D + 1)
    for (i, j), lm in zip(P.monomials, sk.lam):
        e = q ** i + q ** j
        F[e] = K.add(F[e], lm)
    if consts is not None:                       # Dragon: L_k(y) X^(q^k)
        for k in range(d + 1):
            F[q ** k] = K.add(F[q ** k], consts[k])
    for k in range(d):                           # IP: bilinear in X and z
        acc = K.ZERO
        for j in range(r):
            if z[j]:
                acc = K.add(acc, K.scal(sk.Mbil[k][j], z[j]))
        F[q ** k] = K.add(F[q ** k], acc)
    cst = K.ZERO                                 # IP: quadratic in z alone
    for i in range(r):
        for j in range(i, r):
            c = fq.mul(z[i], z[j])
            if c:
                cst = K.add(cst, K.scal(sk.G[i][j], c))
    F[0] = K.add(F[0], cst)
    return F


def _try_solve(P, sk, F, edf):
    """Root-find and check the IP guess; returns the signature or None."""
    K, fq = sk.K, sk.K.fq
    for X in roots(K, norm(K, list(F)), edf):
        if K.is_zero(X):
            continue                             # homogeneous: 0 is trivial
        x = K.coords(X)
        yield x


def sign(P, sk, msg):
    """Sign msg.  Deterministic: no entropy beyond the key and the message."""
    K, fq, q, r, ny = sk.K, sk.K.fq, P.q, P.r, P.ny
    edf = XOF(DOM_EDF, sk.prf, msg)
    zs = list(itertools.product(range(q), repeat=r))

    for counter in range(MAX_SALT):
        if ny is not None:
            # D-James: the salt drives the hash, and the system is solved to 0.
            y = hash_message(msg, counter, ny, q)
            consts = []
            for k in range(P.d + 1):
                acc = K.ZERO
                for j in range(ny):
                    if y[j]:
                        acc = K.add(acc, K.scal(sk.Lam[k][j], y[j]))
                consts.append(acc)
            target = None
        else:
            # James: no salt; the a random minus symbols supply the randomness.
            consts = None
            h = hash_message(msg, 0, P.m, q)
            pad = sample_fq(XOF(DOM_PAD, sk.prf, msg,
                                counter.to_bytes(8, "little")), P.a, q)
            c = list(h) + list(pad)
            u = mat_vec(fq, sk.Tinv, c)
            target = K.from_coords(u)

        for z in zs:
            F = _central_poly(P, sk, consts, z)
            if target is not None:
                F[0] = K.sub(F[0], target)
            for x in _try_solve(P, sk, F, edf):
                if mat_vec(fq, sk.MZ, x) != list(z):
                    continue                     # IP guess was wrong
                a = vec_mat(fq, x, sk.Sinv)
                if any(a):
                    return encode_signature(P, a)
    raise RuntimeError("signing failed after %d attempts" % MAX_SALT)


def encode_signature(P, a):
    """n symbols of F_q, canonically, in exactly ceil(n log2 q / 8) bytes."""
    return encode_vector(a, P.q)


def decode_signature(P, sig):
    """Raises ValueError on a wrong length or a non-canonical encoding."""
    return decode_vector(sig, P.n, P.q)


# ----------------------------------------------------------- verification


def _evaluate(P, pk, a, b):
    """The m public equations at (a, b), as one packed F_q^m vector."""
    W, fq, n, ny = pk.W, None, P.n, P.ny
    from .ff import Fq
    fq = Fq(P.q)
    off = pk._off
    acc = W.ZERO
    k = 0
    nz = [i for i in range(n) if a[i]]
    for i in nz:
        ai = a[i]
        base = off[i]
        for j in nz:
            if j < i:
                continue
            c = fq.mul(ai, a[j])
            if c:
                acc = W.add(acc, W.scal(pk.aa[base + j - i], c))
                k += 1
                if k % _REDUCE_EVERY == 0:
                    acc = W.reduce(acc)
    if ny:
        for i in nz:
            ai, base = a[i], i * ny
            for j in range(ny):
                if b[j]:
                    c = fq.mul(ai, b[j])
                    if c:
                        acc = W.add(acc, W.scal(pk.ab[base + j], c))
                        k += 1
                        if k % _REDUCE_EVERY == 0:
                            acc = W.reduce(acc)
    return W.reduce(acc)


def verify(P, pk, msg, sig):
    try:
        a = decode_signature(P, sig)
    except ValueError:
        return False
    if not any(a):
        return False                             # the trivial root
    W, q = pk.W, P.q
    if P.ny is not None:
        for salt in range(MAX_SALT):
            b = hash_message(msg, salt, P.ny, q)
            if W.is_zero(_evaluate(P, pk, a, b)):
                return True
        return False
    h = hash_message(msg, 0, P.m, q)
    return _evaluate(P, pk, a, None) == W.from_list(h)

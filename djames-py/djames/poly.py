"""Univariate polynomials over K = F_q^n, and root finding.

This is the trapdoor: inverting the D-James central map means finding a root
in K of a polynomial of degree D = q^d + 1, which the parameter sets keep at
24 or below.  Everything here therefore optimises for "tiny degree, huge
coefficient field" -- coefficient lists, schoolbook products, and a Frobenius
step that never touches an exponent larger than q.

Polynomials are plain Python lists of K elements, index = degree, with no
trailing zero coefficient.  The zero polynomial is [].
"""

from .symmetric import XOF, DOM_EDF


def deg(f):
    return len(f) - 1


def norm(K, f):
    while f and K.is_zero(f[-1]):
        f.pop()
    return f


def add(K, a, b):
    n = max(len(a), len(b))
    out = [K.ZERO] * n
    for i, c in enumerate(a):
        out[i] = c
    for i, c in enumerate(b):
        out[i] = K.add(out[i], c)
    return norm(K, out)


def sub(K, a, b):
    n = max(len(a), len(b))
    out = [K.ZERO] * n
    for i, c in enumerate(a):
        out[i] = c
    for i, c in enumerate(b):
        out[i] = K.sub(out[i], c)
    return norm(K, out)


def scal(K, a, s):
    if K.is_zero(s):
        return []
    return norm(K, [K.mul(c, s) for c in a])


def mul(K, a, b):
    if not a or not b:
        return []
    out = [K.ZERO] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if K.is_zero(x):
            continue
        for j, y in enumerate(b):
            if not K.is_zero(y):
                out[i + j] = K.add(out[i + j], K.mul(x, y))
    return norm(K, out)


def monic(K, a):
    if not a:
        return []
    if K.eq(a[-1], K.ONE):
        return a
    inv = K.inv(a[-1])
    return [K.mul(c, inv) for c in a]


def divmod_(K, a, b):
    """(quotient, remainder) of a by b."""
    if not b:
        raise ZeroDivisionError("division by the zero polynomial")
    db = deg(b)
    binv = K.inv(b[-1])
    r = list(a)
    q = [K.ZERO] * max(0, len(a) - db)
    while norm(K, r) and deg(r) >= db:
        dr = deg(r)
        c = K.mul(r[-1], binv)
        q[dr - db] = c
        for j in range(db + 1):
            r[dr - db + j] = K.sub(r[dr - db + j], K.mul(c, b[j]))
        r = norm(K, r)
    return norm(K, q), norm(K, r)


def rem(K, a, b):
    return divmod_(K, a, b)[1]


def gcd(K, a, b):
    a, b = norm(K, list(a)), norm(K, list(b))
    while b:
        a, b = b, rem(K, a, b)
    return monic(K, a)


def powmod(K, base, e, F):
    """base^e mod F."""
    r, b = [K.ONE], rem(K, base, F)
    while e:
        if e & 1:
            r = rem(K, mul(K, r, b), F)
        e >>= 1
        if e:
            b = rem(K, mul(K, b, b), F)
    return r


def evaluate(K, f, x):
    acc = K.ZERO
    for c in reversed(f):
        acc = K.add(K.mul(acc, x), c)
    return acc


# ------------------------------------------------------------ Frobenius on R
#
# The costly part of root finding is X^(q^n) mod F.  We never exponentiate by
# q^n directly; instead we compose the q-power Frobenius with itself, doubling
# the exponent each step, so the cost is O(log n) polynomial products rather
# than O(n).


def _compose_frob(K, F, A, a, B):
    """Given A = X^(q^a) mod F and B = X^(q^b) mod F, return X^(q^(a+b)) mod F.

    B(X)^(q^a) = sum_i B_i^(q^a) * (X^(q^a))^i, and X^(q^a) = A mod F.
    """
    if not B:
        return []
    pw = [[K.ONE]]
    for _ in range(deg(B)):
        pw.append(rem(K, mul(K, pw[-1], A), F))
    out = []
    for i, c in enumerate(B):
        if not K.is_zero(c):
            out = add(K, out, scal(K, pw[i], K.frob(c, a)))
    return out


def x_pow_qn(K, F):
    """X^(q^n) mod F, where n = [K : F_q]."""
    n = K.n
    Q = rem(K, [K.ZERO] * K.q + [K.ONE], F)     # X^q mod F
    res, res_exp = None, 0
    for bit in bin(n)[2:]:
        if res is None:
            res, res_exp = Q, 1                 # leading bit is always 1
        else:
            res = _compose_frob(K, F, res, res_exp, res)
            res_exp *= 2
            if bit == "1":
                res = _compose_frob(K, F, Q, 1, res)
                res_exp += 1
    assert res_exp == n
    return res


# ------------------------------------------------------------- root  finding


def _split(K, G, xof):
    """Factor a squarefree product of distinct linear factors into its roots."""
    if deg(G) == 1:
        return [K.neg(G[0])]                    # G is monic: X + G0
    q, n = K.q, K.n
    while True:
        delta = K.random_nonzero(xof)
        if K.fq.p == 2:
            # char 2: the absolute trace F_{2^m} -> F_2 separates the roots.
            m = n * K.fq.k
            h, term = [], [K.ZERO, delta]       # delta * X
            for _ in range(m):
                h = add(K, h, term)
                term = rem(K, mul(K, term, term), G)
            c = gcd(K, G, h)
        else:
            h = powmod(K, [delta, K.ONE], (q ** n - 1) // 2, G)
            c = gcd(K, G, sub(K, h, [K.ONE]))
        if 0 < deg(c) < deg(G):
            other, r = divmod_(K, G, c)
            assert not r
            return _split(K, monic(K, c), xof) + _split(K, monic(K, other), xof)


def roots(K, F, xof=None):
    """All roots of F in K, as a list (possibly empty).

    gcd(X^(q^n) - X, F) strips F down to the product of its distinct linear
    factors; equal-degree splitting then peels those apart.
    """
    F = monic(K, norm(K, list(F)))
    if deg(F) < 1:
        return []
    if deg(F) == 1:
        return [K.neg(F[0])]
    A = x_pow_qn(K, F)
    G = monic(K, gcd(K, F, sub(K, A, [K.ZERO, K.ONE])))
    if deg(G) < 1:
        return []
    if xof is None:
        xof = XOF(DOM_EDF)
    return _split(K, G, xof)

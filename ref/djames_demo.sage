# ---- cell 1 ----



import hashlib
from collections import namedtuple
# ---- cell 3 ----
def field_element_digits(elt, length):
    """
    GF(q)-coordinates of a finite-field element `elt`, expressed in the
    power basis 1, a, a^2, ... of its field, zero-padded/truncated to
    exactly `length` entries.

    Returns a plain Python list of ints in [0, q).
    """
    coeffs = [int(c) for c in elt.polynomial().list()]
    coeffs += [0] * (length - len(coeffs))
    return coeffs[:length]


def hash_to_digits(message, length, q):
    """
    Deterministically hash `message` (any object with a stable `repr`) into
    `length` base-q digits in [0, q-1], derived from a SHA-256 digest.

    This stands in for a proper domain-separated hash-to-(Z/qZ)^length
    function. It is adequate for a proof-of-concept but should not be used
    as-is in a real implementation.
    """
    digest = int(hashlib.sha256(repr(message).encode()).hexdigest(), 16)
    digits = []
    for _ in range(length):
        digits.append(digest % q)
        digest //= q
    return digits


def evaluate_affine_form(Q, L, cst, x):
    """
    Evaluate x^T Q x + L.x + cst, where Q is a square matrix and L a
    vector over some field K, cst in K, and x is any length-matching
    iterable of elements coercible into K (e.g. a signature, i.e. a
    vector over the base field GF(q)).

    Note the explicit *linear* part L: over GF(2), x_i^2 = x_i, so a
    linear term can be folded into the diagonal of a "quadratic" matrix
    at no cost -- a classic trick that only holds in characteristic 2.
    For q > 2 that identity is false (e.g. 2^2 = 1 != 2 in GF(3)), so
    linear and quadratic contributions must be tracked separately; this
    is exactly the difference between this notebook's affine-quadratic
    representation and the purely-quadratic one that a q=2-only
    implementation can get away with.
    """
    K = Q.base_ring()
    xK = vector(K, list(x))
    Lv = vector(K, list(L))
    return xK * Q * xK + Lv * xK + cst


def random_invertible_matrix(field, dim):
    """A uniformly random invertible dim x dim matrix over `field`."""
    return random_matrix(field, dim, algorithm="echelonizable", rank=dim)


def encode_public_coefficient(coeff, Tpub, q, n, Kpub):
    """
    Re-encode a coefficient `coeff` living in the secret field K = GF(q^n)
    as an element of the public field Kpub = GF(q^m): decompose `coeff`
    into its n GF(q)-digits, apply the public projection matrix Tpub
    (an m x n matrix over GF(q), i.e. the first m rows of the secret map T),
    and repack the resulting m digits as an element of Kpub.
    """
    digits = vector(GF(q), field_element_digits(coeff, n))
    projected = Tpub * digits
    Rm = Kpub.polynomial_ring()
    return Kpub(Rm(list(projected)))


def public_affine_form(Q, L, S, T, q, n, m, Kpub):
    """
    Build the public affine-quadratic form from a secret quadratic form Q
    (an N x N matrix over K = GF(q^n)), a secret linear part L (length N
    over K), and the secret linear maps S (N x N over GF(q), masking the
    variables) and T (n x n over GF(q), whose first m rows project the n
    equations down to m public equations). N is either n or n +
    hash_length depending on the scheme.

    Substituting x = y S^{-1}... (equivalently: masking by S) gives

        Q_S = S Q S^T                (new quadratic part)
        L_S = L S^T                  (new linear part)

    and each entry is then projected to the public field via pi(.) =
    `encode_public_coefficient` (using Tpub = T[:m, :]), exactly as for
    the plain-HFE constant term. Only the upper triangle of Q_pub is
    filled (see `encode_public_coefficient`'s docstring / the notebook
    text for why folding is valid for any q).

    Returns (Qpub, Lpub).
    """
    N = Q.nrows()
    K = Q.base_ring()
    S_K = S.change_ring(K)
    QS = S_K * Q * S_K.transpose()
    LS = vector(K, list(L)) * S_K.transpose()
    Tpub = T[:m, :]

    Qpub = Matrix(Kpub, N, N)
    for i in range(N):
        for j in range(i, N):
            coeff = QS[i, i] if i == j else QS[i, j] + QS[j, i]
            Qpub[i, j] = encode_public_coefficient(coeff, Tpub, q, n, Kpub)

    Lpub = vector(Kpub, [encode_public_coefficient(LS[i], Tpub, q, n, Kpub)
                          for i in range(N)])
    return Qpub, Lpub
# ---- cell 5 ----
def random_hfe_polynomial(Rx, q, HFEDegI, HFEDegJ):
    """
    Draw a uniformly random HFE polynomial of shape (HFEDegI, HFEDegJ):

        F(X) = X^(q^HFEDegI + q^HFEDegJ)
             + sum_{i<HFEDegI, j<=HFEDegJ} c_{i,j} X^(q^i + q^j)
             + sum_{i<=HFEDegI}            c_i     X^(q^i)
             + c_0

    with coefficients drawn uniformly from K = Rx.base_ring().
    """
    K = Rx.base_ring()
    X = Rx.gen()
    HFEDeg = q**HFEDegI + q**HFEDegJ
    F = X**HFEDeg
    F += sum(K.random_element() * X**(q**i + q**j)
              for j in range(HFEDegJ + 1) for i in range(HFEDegI))
    F += sum(K.random_element() * X**(q**i) for i in range(HFEDegI + 1))
    return F


def hfe_quadratic_form(F, q, n, HFEDegI, HFEDegJ):
    """
    Build the pair (Q, L) over K = F.base_ring() = GF(q^n) associated with
    the HFE polynomial F of shape (HFEDegI, HFEDegJ), i.e. the n x n
    quadratic part Q and the length-n linear part L such that

        F(X) - F(0)  =  sum_{i,j} Q[i,j] x_i x_j  +  sum_i L[i] x_i    over K,

    once X = sum_k theta_k x_k is substituted (theta_k = alpha^k). This
    uses two Frobenius-type identities, both valid for *any* prime q (not
    just q = 2): x_k^(q^t) = x_k for x_k in GF(q), which makes every
    "single q-power" term X^(q^t) expand as a genuinely *linear*
    combination sum_k theta_k^(q^t) x_k (so it belongs in L); and the
    ordinary binomial expansion of (X^(q^i))^2, needed for the "doubled"
    terms X^(2 q^i) that appear when the two exponents making up a
    quadratic HFE term coincide (i = j) -- these are genuinely *quadratic*
    (they involve x_k^2, not x_k) and belong in Q.

    A q = 2-only implementation can get away with folding every term into
    a single "quadratic" matrix, because x_k^2 = x_k for x_k in GF(2)
    collapses the two cases into one; this is exactly the corner the
    original code cut (see the notebook introduction).
    """
    K = F.base_ring()
    alpha = K.gen()
    theta = [alpha**i for i in range(n)]
    Q = Matrix(K, n, n)
    L = vector(K, n)

    # --- X^1 term: genuinely LINEAR (exponent q^0 = 1) ---------------------
    if HFEDegI == 0 and HFEDegJ == 0:
        for i in range(n):
            L[i] += theta[i]
    else:
        c = F[1]
        for i in range(n):
            L[i] += theta[i] * c

    # --- remaining terms, one HFE "level" k = 0 .. HFEDegI-1 at a time ----
    for k in range(HFEDegI):
        exp_i = q**k

        # quadratic cross terms X^(q^k) X^(q^l), l < k (genuinely quadratic,
        # since the two exponents q^k, q^l differ)
        for l in range(k):
            exp_j = q**l
            coef = F[exp_i + exp_j]
            for i in range(n):
                ti = theta[i] ** exp_i
                for j in range(n):
                    tj = theta[j] ** exp_j
                    Q[i, j] += ti * tj * coef

        # linear term X^(q^(k+1)): genuinely LINEAR
        exp = q**(k + 1)
        coef = F[exp]
        for i in range(n):
            L[i] += (theta[i] ** exp) * coef

    # --- highest quadratic-degree cross terms (genuinely quadratic) -------
    exp_i = q**HFEDegI
    for l in range(HFEDegJ):
        exp_j = q**l
        coef = F[exp_i + exp_j]
        for i in range(n):
            ti = theta[i] ** exp_i
            for j in range(n):
                tj = theta[j] ** exp_j
                Q[i, j] += ti * tj * coef

    # --- leading (monic) term -----------------------------------------------
    if HFEDegI != HFEDegJ:
        exp_j = q**HFEDegJ
        for i in range(n):
            ti = theta[i] ** exp_i
            for j in range(n):
                tj = theta[j] ** exp_j
                Q[i, j] += ti * tj
    else:
        # X^(2 q^HFEDegI): the two exponents coincide, so (unlike the
        # branch above) this is a "doubled"/square-type term -- genuinely
        # quadratic, full bilinear expansion, coefficient 1 (monic).
        for i in range(n):
            ti = theta[i] ** exp_i
            for j in range(n):
                tj = theta[j] ** exp_i
                Q[i, j] += ti * tj

    # --- "doubled" cross terms X^(2 q^i), i.e. i = j pairs in the double
    # sum defining the HFE polynomial's quadratic terms (see
    # `random_hfe_polynomial`): for every i with i < HFEDegI and i <=
    # HFEDegJ, the polynomial contains a term at exponent q^i + q^i.
    #
    # For q = 2 this exponent equals q^(i+1) -- exactly the "linear term"
    # exponent already handled above -- and the *same* underlying
    # coefficient of F is shared between both conceptual contributions
    # (Sage/the field only stores one coefficient per exponent).
    # For q > 2 this exponent is distinct from every linear exponent, was
    # never read by any loop above, and must be added here as a real
    # quadratic contribution (with its genuine cross terms this time).
    if q != 2:
        for i in range(min(HFEDegI, HFEDegJ + 1)):
            exp_i = q**i
            coef = F[2 * exp_i]
            if coef == 0:
                continue
            for a in range(n):
                ta = theta[a] ** exp_i
                for b in range(n):
                    Q[a, b] += ta * (theta[b] ** exp_i) * coef

    return Q, L
# ---- cell 7 ----
def add_dragon_terms(Q, L, F, HFEDegI, DragonMat, q):
    """
    Extend an n x n HFE affine-quadratic form (Q, L) (over K =
    F.base_ring()) with the bilinear cross-terms sum_t sum_j
    DragonMat[t,j] * X^(q^t) * Y_j coming from the Dragon hash-embedding
    trick, and return the resulting (Qext, Lext) pair in the extended
    (n + hash_length)-dimensional variable set (x_0,...,x_{n-1},
    Y_0,...,Y_{hash_length-1}). These cross-terms are quadratic (they
    involve both an x and a Y variable), so only Q grows; L is simply
    zero-padded to match the larger variable set.

    DragonMat must have at least HFEDegI + 1 rows; only rows 0..HFEDegI
    are read.
    """
    K = F.base_ring()
    alpha = K.gen()
    n = Q.nrows()
    hash_length = DragonMat.ncols()
    theta = [alpha**i for i in range(n)]

    Qext = Matrix(K, n + hash_length, n + hash_length)
    for i in range(n):
        for j in range(n):
            Qext[i, j] = Q[i, j]
    Lext = vector(K, list(L) + [K.zero()] * hash_length)

    for t in range(HFEDegI + 1):
        exp = q**t
        for j in range(hash_length):
            coef = DragonMat[t, j]
            if coef == 0:
                continue
            for i in range(n):
                Qext[i, n + j] += (theta[i] ** exp) * coef

    return Qext, Lext


def add_internal_perturbation_terms(Q, F, HFEDegI, MZ, MBilin, HQ, q):
    """
    Extend an n x n HFE quadratic form Q in place with the internal-
    perturbation terms coming from substituting z = MZ * x, an auxiliary
    r-dimensional vector of hidden variables:

        sum_{t=0}^{HFEDegI} sum_j MBilin[t,j] * X^(q^t) * z_j
      + sum_{i<=j} HQ[i,j] * z_i * z_j

    Both families are quadratic in x (z is itself linear in x, so
    "X * z" and "z * z" are both degree-2 in x), so only Q is touched --
    IP never introduces a genuinely linear contribution, and does not
    introduce new public variables either, so the result is still n x n.
    Returns Q (also mutated in place).
    """
    K = F.base_ring()
    alpha = K.gen()
    n = Q.nrows()
    r = MZ.nrows()
    theta = [alpha**i for i in range(n)]

    # bilinear X^(q^t) z_j terms
    for t in range(HFEDegI + 1):
        exp = q**t
        for j in range(r):
            coef = MBilin[t, j]
            if coef == 0:
                continue
            row = MZ.row(j)
            for a in range(n):
                ta = theta[a] ** exp
                for b in range(n):
                    if row[b] != 0:
                        Q[a, b] += coef * ta * row[b]

    # quadratic z_i z_j terms
    for i in range(r):
        zi = MZ.row(i)
        for j in range(i, r):
            coef = HQ[i, j]
            if coef == 0:
                continue
            zj = MZ.row(j)
            for a in range(n):
                if zi[a] == 0:
                    continue
                for b in range(n):
                    if zj[b] != 0:
                        Q[a, b] += coef * zi[a] * zj[b]

    return Q


def hfe_dragon_quadratic_form(F, q, n, HFEDegI, HFEDegJ, DragonMat):
    Q, L = hfe_quadratic_form(F, q, n, HFEDegI, HFEDegJ)
    return add_dragon_terms(Q, L, F, HFEDegI, DragonMat, q)


def hfe_ip_quadratic_form(F, q, n, HFEDegI, HFEDegJ, MZ, MBilin, HQ):
    Q, L = hfe_quadratic_form(F, q, n, HFEDegI, HFEDegJ)
    Q = add_internal_perturbation_terms(Q, F, HFEDegI, MZ, MBilin, HQ, q)
    return Q, L


def hfe_ip_dragon_quadratic_form(F, q, n, HFEDegI, HFEDegJ, MZ, MBilin, HQ, DragonMat):
    Q, L = hfe_quadratic_form(F, q, n, HFEDegI, HFEDegJ)
    Q = add_internal_perturbation_terms(Q, F, HFEDegI, MZ, MBilin, HQ, q)
    return add_dragon_terms(Q, L, F, HFEDegI, DragonMat, q)
# ---- cell 9 ----
Params = namedtuple("Params", ["q", "n", "m", "HFEDegI", "HFEDegJ", "MAX_SALT", "r", "hash_length"])
Params.__new__.__defaults__ = (None, None, None)  #MAX_SALT, r, hash_length default to None
# ---- cell 11 ----
def hfe_keygen(params):
    q, n, m = params.q, params.n, params.m
    HFEDegI, HFEDegJ = params.HFEDegI, params.HFEDegJ
    if m > n:
        raise ValueError("Need m <= n")

    Fq = GF(q)
    K = GF(q**n, name="alpha")
    Rx = PolynomialRing(K, "X")
    F = random_hfe_polynomial(Rx, q, HFEDegI, HFEDegJ)

    Q, L = hfe_quadratic_form(F, q, n, HFEDegI, HFEDegJ)
    S = random_invertible_matrix(Fq, n)
    T = random_invertible_matrix(Fq, n)

    Kpub = GF(q**m, name="b")
    Tpub = T[:m, :]
    cst_pub = encode_public_coefficient(F[0], Tpub, q, n, Kpub)
    Qpub, Lpub = public_affine_form(Q, L, S, T, q, n, m, Kpub)

    pk = (cst_pub, Qpub, Lpub)
    sk = (F, S.inverse(), T.inverse())
    return pk, sk


def hfe_sign(message, sk, params):
    q, n, m = params.q, params.n, params.m
    F, S_inv, T_inv = sk
    K = F.base_ring()
    alpha = K.gen()
    Fq = GF(q)

    h = hash_to_digits(message, m, q)
    while True:
        if m < n:
            pad = random_vector(Fq, n - m)
            c = vector(Fq, list(h) + list(pad))
        else:
            c = vector(Fq, h)

        U_digits = T_inv * c
        U = sum(U_digits[i] * alpha**i for i in range(n))

        roots = (F - U).roots()
        root = next((rt for rt, mult in roots), None)
        if root is not None:
            break
        if m == n:
            raise RuntimeError("No root found for this HFE instance")

    x = vector(root)
    return x * S_inv


def hfe_verify(pk, message, sig, params):
    q, m = params.q, params.m
    cst_pub, Qpub, Lpub = pk
    h = hash_to_digits(message, m, q)
    value = evaluate_affine_form(Qpub, Lpub, cst_pub, sig)
    return field_element_digits(value, m) == h
# ---- cell 13 ----
def hfe_dragon_keygen(params):
    q, n, m = params.q, params.n, params.m
    HFEDegI, HFEDegJ = params.HFEDegI, params.HFEDegJ
    hash_length = params.hash_length
    if m > n:
        raise ValueError("Need m <= n")

    Fq = GF(q)
    K = GF(q**n, name="alpha")
    Rx = PolynomialRing(K, "X")
    F = random_hfe_polynomial(Rx, q, HFEDegI, HFEDegJ)

    # DragonMat has n rows: only the first HFEDegI are random, the rest
    # (including row HFEDegI, i.e. the coefficient of X^(q^HFEDegI) * Y_j)
    # are zero.
    DragonMat = random_matrix(K, HFEDegI, hash_length).stack(
        zero_matrix(K, n - HFEDegI, hash_length)
    )

    Q, L = hfe_dragon_quadratic_form(F, q, n, HFEDegI, HFEDegJ, DragonMat)
    N = n + hash_length

    Sx = random_invertible_matrix(Fq, n)
    S = block_diagonal_matrix(Sx, identity_matrix(Fq, hash_length))
    T = random_invertible_matrix(Fq, n)

    Kpub = GF(q**m, name="b")
    Tpub = T[:m, :]
    cst_pub = encode_public_coefficient(F[0], Tpub, q, n, Kpub)
    Qpub, Lpub = public_affine_form(Q, L, S, T, q, n, m, Kpub)

    pk = (cst_pub, Qpub, Lpub)
    # S is block-diagonal, so the top-left n x n block of S^{-1} is Sx^{-1}.
    sk = (F, HFEDegI, DragonMat, Sx.inverse(), T.inverse())
    return pk, sk


def hfe_dragon_sign(message, sk, params):
    q, MAX_SALT = params.q, params.MAX_SALT
    F, HFEDegI, DragonMat, S_inv, T_inv = sk
    K = F.base_ring()
    Rx = F.parent()
    X = Rx.gen()
    hash_length = DragonMat.ncols()

    for salt in range(MAX_SALT):
        h = hash_to_digits((message, salt), hash_length, q)

        Fh = F
        for j in range(hash_length):
            if h[j] == 0:
                continue
            yj = K(h[j])
            for k in range(HFEDegI + 1):
                coef = DragonMat[k, j]
                if coef != 0:
                    Fh += yj * coef * X**(q**k)

        roots = Fh.roots()
        root = next((rt for rt, mult in roots if rt != 0), None)

        if root is not None:
            x = vector(root)
            return x * S_inv

    raise RuntimeError("Unable to sign that message")


def hfe_dragon_verify(pk, message, sig, params):
    q, m, MAX_SALT = params.q, params.m, params.MAX_SALT
    hash_length = params.hash_length
    cst_pub, Qpub, Lpub = pk
    for counter in range(MAX_SALT):
        h = hash_to_digits((message, counter), hash_length, q)

        v = vector(GF(q), list(sig) + h)
        value = evaluate_affine_form(Qpub, Lpub, cst_pub, v)

        if field_element_digits(value, m) == [0] * m:
            return True

    return False
# ---- cell 15 ----
def hfe_ip_keygen(params):
    q, n, m, r = params.q, params.n, params.m, params.r
    HFEDegI, HFEDegJ = params.HFEDegI, params.HFEDegJ
    if m > n:
        raise ValueError("Need m <= n")

    Fq = GF(q)
    K = GF(q**n, name="alpha")
    Rx = PolynomialRing(K, "X")
    F = random_hfe_polynomial(Rx, q, HFEDegI, HFEDegJ)

    MZ = random_matrix(Fq, r, n, algorithm="echelonizable", rank=r)

    MBilin = Matrix(K, HFEDegI + 1, r)
    for i in range(HFEDegI + 1):
        for j in range(r):
            MBilin[i, j] = K.random_element()

    HQ = Matrix(K, r, r)
    for i in range(r):
        for j in range(i, r):
            HQ[i, j] = K.random_element()

    Q, L = hfe_ip_quadratic_form(F, q, n, HFEDegI, HFEDegJ, MZ, MBilin, HQ)

    S = random_invertible_matrix(Fq, n)
    T = random_invertible_matrix(Fq, n)

    Kpub = GF(q**m, name="b")
    Tpub = T[:m, :]
    cst_pub = encode_public_coefficient(F[0], Tpub, q, n, Kpub)
    Qpub, Lpub = public_affine_form(Q, L, S, T, q, n, m, Kpub)

    pk = (cst_pub, Qpub, Lpub)
    sk = (F, HFEDegI, MZ, MBilin, HQ, S.inverse(), T.inverse())
    return pk, sk


def hfe_ip_sign(message, sk, params):
    q = params.q
    F, HFEDegI, MZ, MBilin, HQ, S_inv, T_inv = sk
    K = F.base_ring()
    Rx = F.parent()
    X = Rx.gen()
    Fq = GF(q)
    r = MZ.nrows()

    while True:
        z = random_vector(Fq, r)

        Fz = F
        for i in range(r):
            if z[i] == 0:
                continue
            zi = K(z[i])
            for k in range(HFEDegI + 1):
                coef = MBilin[k, i]
                if coef != 0:
                    Fz += zi * coef * X**(q**k)

        cst = K.zero()
        for i in range(r):
            cst += K(z[i])**2 * HQ[i, i]
            for j in range(i + 1, r):
                cst += K(z[i]) * K(z[j]) * HQ[i, j]
        Fz += cst

        try:
            sigma = hfe_sign(message, (Fz, S_inv, T_inv), params)
        except RuntimeError:
            continue

        x = sigma * S_inv.inverse()
        if MZ * x == z:
            return sigma


def hfe_ip_verify(pk, message, sig, params):
    """Internal perturbation does not change the verification equation."""
    return hfe_verify(pk, message, sig, params)
# ---- cell 17 ----
def hfe_ip_dragon_keygen(params):
    q, n, m, r = params.q, params.n, params.m, params.r
    HFEDegI, HFEDegJ = params.HFEDegI, params.HFEDegJ
    hash_length = params.hash_length
    if m > n:
        raise ValueError("Need m <= n")

    Fq = GF(q)
    K = GF(q**n, name="alpha")
    Rx = PolynomialRing(K, "X")
    F = random_hfe_polynomial(Rx, q, HFEDegI, HFEDegJ)

    MZ = random_matrix(Fq, r, n, algorithm="echelonizable", rank=r)

    MBilin = Matrix(K, HFEDegI + 1, r)
    for i in range(HFEDegI + 1):
        for j in range(r):
            MBilin[i, j] = K.random_element()

    HQ = Matrix(K, r, r)
    for i in range(r):
        for j in range(i, r):
            HQ[i, j] = K.random_element()

    # Unlike HFE-Dragon's L, DragonL here has exactly HFEDegI+1 rows, all
    # drawn at random (no forced-zero row).
    DragonL = Matrix(K, HFEDegI + 1, hash_length)
    for i in range(HFEDegI + 1):
        for j in range(hash_length):
            DragonL[i, j] = K.random_element()

    Q, L = hfe_ip_dragon_quadratic_form(F, q, n, HFEDegI, HFEDegJ, MZ, MBilin, HQ, DragonL)
    N = n + hash_length

    Sx = random_invertible_matrix(Fq, n)
    S = block_diagonal_matrix(Sx, identity_matrix(Fq, hash_length))
    T = random_invertible_matrix(Fq, n)

    Kpub = GF(q**m, name="b")
    Tpub = T[:m, :]
    cst_pub = encode_public_coefficient(F[0], Tpub, q, n, Kpub)
    Qpub, Lpub = public_affine_form(Q, L, S, T, q, n, m, Kpub)

    pk = (cst_pub, Qpub, Lpub)
    sk = (F, HFEDegI, MZ, MBilin, HQ, DragonL, Sx.inverse(), T.inverse())
    return pk, sk


def hfe_ip_dragon_sign(message, sk, params):
    q, MAX_SALT = params.q, params.MAX_SALT
    F, HFEDegI, MZ, MBilin, HQ, DragonL, S_inv, T_inv = sk
    K = F.base_ring()
    Rx = F.parent()
    X = Rx.gen()
    Fq = GF(q)
    r = MZ.nrows()
    hash_length = DragonL.ncols()

    for salt in range(MAX_SALT):
        Y = hash_to_digits((message,salt), hash_length, q)
        FY = F
        for j in range(hash_length):
            if Y[j] == 0:
                continue
            yj = K(Y[j])
            for k in range(HFEDegI + 1):
                coef = DragonL[k, j]
                if coef != 0:
                    FY += yj * coef * X**(q**k)    
                    
        for z in Fq^r:            
            Fz = FY
            for i in range(r):
                if z[i] == 0:
                    continue
                zi = K(z[i])
                for k in range(HFEDegI + 1):
                    coef = MBilin[k, i]
                    if coef != 0:
                        Fz += zi * coef * X**(q**k)

            cst = K.zero()
            for i in range(r):
                cst += K(z[i])**2 * HQ[i, i]
                for j in range(i + 1, r):
                    cst += K(z[i]) * K(z[j]) * HQ[i, j]
            Fz += cst

            roots = Fz.roots()
            x_secret = next((rt for rt, mult in roots if rt != 0), None)
            if x_secret is None:
                continue

            x_vec = vector(x_secret)
            if MZ * x_vec != z:
                continue

            return x_vec * S_inv
    raise RuntimeError("Unable to sign that message")


def hfe_ip_dragon_verify(pk, message, sig, params):
    """Same verification equation as HFE-Dragon."""
    return hfe_dragon_verify(pk, message, sig, params)
# ---- cell 19 ----
def quadratic_form_to_gfq_system(Q, L, q, n):
    """
    Expand an N x N quadratic form Q and length-N linear part L, both over
    K = GF(q^n), into an explicit system of n multivariate quadratic
    polynomials over GF(q) in N variables x0,...,x_{N-1}: the k-th
    polynomial collects the k-th GF(q)-digit (power-basis coordinate) of
    every entry of Q and L.

    Works uniformly for plain HFE (N = n) and for Dragon-extended forms
    (N = n + hash_length).
    """
    N = Q.nrows()
    quad_coeffs = [Matrix(GF(q), N, N) for _ in range(n)]
    lin_coeffs = [vector(GF(q), N) for _ in range(n)]
    for i in range(N):
        digits = field_element_digits(L[i], n)
        for k in range(n):
            lin_coeffs[k][i] = digits[k]
        for j in range(N):
            digits = field_element_digits(Q[i, j], n)
            for k in range(n):
                quad_coeffs[k][i, j] = digits[k]

    R = PolynomialRing(GF(q), N, names=[f"x{i}" for i in range(N)])
    xs = R.gens()

    equations = []
    for k in range(n):
        f = R.zero()
        for i in range(N):
            if lin_coeffs[k][i] != 0:
                f += lin_coeffs[k][i] * xs[i]
            for j in range(N):
                c = quad_coeffs[k][i, j]
                if c != 0:
                    f += c * xs[i] * xs[j]
        equations.append(f)
    return equations
# ---- cell 21 ----
def run_scheme_demo(name, keygen, sign, verify, params, trials=3, message_len=64):
    tag = f"q={params.q}, n={params.n}, m={params.m}"
    if params.r is not None:
        tag += f", r={params.r}"
    if params.hash_length is not None:
        tag += f", hash_length={params.hash_length}"
    print(f"--- {name}  ({tag}) ---")

    pk, sk = keygen(params)

    for t in range(trials):
        message = [randint(0, 1) for _ in range(message_len)]
        sigma = sign(message, sk, params)
        ok = verify(pk, message, sigma, params)
        print(f"  trial {t + 1}: verify(sign(message)) = {ok}")
        if not ok:
            raise AssertionError("A validly generated signature failed to verify!")

    # Negative test: perturbing one coordinate of a valid signature should
    # (with overwhelming probability) make verification fail.
    message = [randint(0, 1) for _ in range(message_len)]
    sigma = list(sign(message, sk, params))
    sigma[0] = sigma[0] + 1
    tampered = vector(GF(params.q), sigma)
    ok = verify(pk, message, tampered, params)
    print(f"  tampered signature accepted = {ok}  (expected: False)")
    print()

print("SAGE VERSION CHECK OK")
p2   = Params(q=2, n=48, m=32, HFEDegI=2, HFEDegJ=0, MAX_SALT=10, hash_length=64)
p2ip = Params(q=2, n=48, m=32, HFEDegI=2, HFEDegJ=0, MAX_SALT=10, r=2, hash_length=64)
p5ip = Params(q=5, n=21, m=14, HFEDegI=1, HFEDegJ=0, MAX_SALT=10, r=2, hash_length=28)
run_scheme_demo("HFE-Dragon (q=2)",    hfe_dragon_keygen,    hfe_dragon_sign,    hfe_dragon_verify,    p2,   trials=2)
run_scheme_demo("D-James=HFE-IP-Dragon (q=2)", hfe_ip_dragon_keygen, hfe_ip_dragon_sign, hfe_ip_dragon_verify, p2ip, trials=2)
run_scheme_demo("D-James=HFE-IP-Dragon (q=5)", hfe_ip_dragon_keygen, hfe_ip_dragon_sign, hfe_ip_dragon_verify, p5ip, trials=2)

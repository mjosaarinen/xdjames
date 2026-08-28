"""Finite fields for D-James: F_q, and the extension K = F_q^n.

Representation
--------------
A straightforward list-of-ints encoding of K is hopeless here: the paper's
parameters reach n = 390 with n_y = 512, and public-key generation touches
O(n * (n + n_y)) field products.  So every vector over F_q -- an element of K,
a matrix row, a public-key coefficient column -- is packed inside a Python
big integer, and the bulk arithmetic is expressed as `int` operations, which
CPython runs at C speed.

Three backends, chosen by q:

  F2Vec    q = 2       one coefficient per bit; addition is XOR
  F2kVec   q = 2^k     k bit-planes; addition is XOR plane-wise
  FpVec    q = p odd   one coefficient per B-bit lane

For odd p the lanes are wide enough that a whole polynomial product is a
*single* integer multiplication: lane k of a*b accumulates sum_{i+j=k} a_i b_j,
which is at most n(q-1)^2, so as long as B exceeds that bit length no lane
can carry into its neighbour.  Reducing the lanes mod p afterwards is the only
per-coefficient work.
"""

# --------------------------------------------------------------- base field


def _is_prime(x):
    if x < 2:
        return False
    d = 2
    while d * d <= x:
        if x % d == 0:
            return False
        d += 1
    return True


def _prime_power(q):
    """Write q = p^k with p prime."""
    for p in range(2, q + 1):
        if q % p:
            continue
        if not _is_prime(p):
            continue
        k, t = 0, q
        while t % p == 0:
            t //= p
            k += 1
        if t != 1:
            break
        return p, k
    raise ValueError("q = %d is not a prime power" % q)


def _first_irreducible_over_fp(p, k):
    """Lex-least monic irreducible of degree k over F_p, as [c0, ..., c_{k-1}].

    Only ever called for tiny (p, k) -- q = 4 is the sole prime power among the
    paper's parameter sets -- so exhaustive trial division is plenty.
    """
    def mul(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] = (r[i + j] + x * y) % p
        return r

    def divides(g, f):
        """Does monic g divide monic f?"""
        f = f[:]
        dg = len(g) - 1
        for i in range(len(f) - 1, dg - 1, -1):
            c = f[i]
            if c:
                for j in range(dg + 1):
                    f[i - dg + j] = (f[i - dg + j] - c * g[j]) % p
        return not any(f[:dg])

    def monics(d):
        for v in range(p ** d):
            cs, t = [], v
            for _ in range(d):
                cs.append(t % p)
                t //= p
            yield cs + [1]

    for cand in monics(k):
        if all(not divides(g, cand) for d in range(1, k // 2 + 1)
               for g in monics(d)):
            return cand[:k]
    raise AssertionError("no irreducible of degree %d over F_%d" % (k, p))


class Fq:
    """The base field F_q.  Elements are ints in [0, q).

    For q = p^k with k > 1 an element encodes its base-p digit vector as the
    coefficients of a polynomial in u, modulo a fixed irreducible.  Full q x q
    operation tables are precomputed; q <= 23 throughout, so this is 529
    entries at worst.
    """

    def __init__(self, q):
        self.q = q
        self.p, self.k = _prime_power(q)
        if self.k == 1:
            self.MUL = [[(a * b) % q for b in range(q)] for a in range(q)]
            self.ADD = [[(a + b) % q for b in range(q)] for a in range(q)]
        else:
            p, k = self.p, self.k
            self.modulus = _first_irreducible_over_fp(p, k)

            def digits(e):
                d, t = [], e
                for _ in range(k):
                    d.append(t % p)
                    t //= p
                return d

            def pack(d):
                return sum(c * p ** i for i, c in enumerate(d))

            def fmul(a, b):
                x, y = digits(a), digits(b)
                r = [0] * (2 * k - 1)
                for i, xi in enumerate(x):
                    if xi:
                        for j, yj in enumerate(y):
                            r[i + j] = (r[i + j] + xi * yj) % p
                for i in range(2 * k - 2, k - 1, -1):   # reduce mod `modulus`
                    c = r[i]
                    if c:
                        r[i] = 0
                        for j in range(k):
                            r[i - k + j] = (r[i - k + j] - c * self.modulus[j]) % p
                return pack(r[:k])

            self.MUL = [[fmul(a, b) for b in range(q)] for a in range(q)]
            self.ADD = [[pack([(x + y) % p for x, y in
                               zip(digits(a), digits(b))]) for b in range(q)]
                        for a in range(q)]

        self.NEG = [self.ADD[a].index(0) for a in range(q)]
        self.INV = [0] * q
        for a in range(1, q):
            self.INV[a] = self.MUL[a].index(1)

    def mul(self, a, b):
        return self.MUL[a][b]

    def add(self, a, b):
        return self.ADD[a][b]

    def sub(self, a, b):
        return self.ADD[a][self.NEG[b]]

    def neg(self, a):
        return self.NEG[a]

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("0 has no inverse in F_q")
        return self.INV[a]


# --------------------------------------------------- packed vectors over F_q
#
# Every backend exposes the same interface.  Lengths are implicit: a packed
# value simply has zero coefficients above its degree, so the same routines
# serve elements of K (length n) and the wide products (length 2n-1) that
# reduction consumes.

_SPREAD16 = [bytes(x for i in range(8) for x in ((b >> i) & 1, 0))
             for b in range(256)]
_PACK8 = {}
for _b in range(256):
    _PACK8[bytes((_b >> _i) & 1 for _i in range(8))] = _b


def _clmul(a, b):
    """Carry-less (F_2[t]) product of two bit-packed polynomials.

    Both operands are spread to 16-bit lanes so that an ordinary integer
    multiply computes every coefficient sum_{i+j=k} a_i b_j without carrying
    between lanes (the sums are bounded by the operand length, far below
    2^16).  Reading the low bit of each lane then recovers the product mod 2.
    """
    if a == 0 or b == 0:
        return 0
    na, nb = (a.bit_length() + 7) // 8, (b.bit_length() + 7) // 8
    A = int.from_bytes(b"".join([_SPREAD16[c] for c in a.to_bytes(na, "little")]), "little")
    B = int.from_bytes(b"".join([_SPREAD16[c] for c in b.to_bytes(nb, "little")]), "little")
    P = A * B
    nl = na + nb                        # bytes of coefficients, i.e. lanes/8
    raw = P.to_bytes(nl * 16, "little")[0::2]     # one byte per coefficient
    raw = bytes(c & 1 for c in raw)
    return int.from_bytes(bytes([_PACK8[raw[i:i + 8]] for i in range(0, len(raw), 8)]),
                          "little")


class F2Vec:
    """Vectors over F_2: coefficient i is bit i of a Python int."""

    q = 2

    def __init__(self):
        self.ZERO = 0

    def add(self, a, b):
        return a ^ b

    sub = add

    def neg(self, a):
        return a

    def scal(self, a, s):
        return a if s else 0

    def lincomb(self, cs, vs):
        r = 0
        for i, c in enumerate(cs):
            if c:
                r ^= vs[i]
        return r

    def reduce(self, a):
        return a

    def polymul(self, a, b):
        return _clmul(a, b)

    def shift(self, a, k):
        return a << k

    def rshift(self, a, k):
        return a >> k

    def truncate(self, a, L):
        return a & ((1 << L) - 1)

    def coef(self, a, i):
        return (a >> i) & 1

    def to_list(self, a, L):
        return [(a >> i) & 1 for i in range(L)]

    def from_list(self, cs):
        r = 0
        for i, c in enumerate(cs):
            if c:
                r |= 1 << i
        return r

    def is_zero(self, a):
        return a == 0

    def degree(self, a):
        return a.bit_length() - 1

    def support(self, a):
        """Indices of the nonzero coefficients -- the fast path for lincomb."""
        out = []
        while a:
            lb = a & -a
            out.append(lb.bit_length() - 1)
            a ^= lb
        return out


class F2kVec:
    """Vectors over F_{2^k}, k > 1: k bit-planes held as a k-tuple of ints.

    Plane j carries the u^j digit of every coefficient, so addition is a
    plane-wise XOR and a product decomposes into k^2 F_2 carry-less products.
    """

    def __init__(self, fq):
        self.fq = fq
        self.q = fq.q
        self.k = fq.k
        self.ZERO = (0,) * self.k
        # mix[s][j] = bit-mask of planes that s * u^j contributes to
        self.mix = [[[(fq.mul(s, 1 << j) >> i) & 1 for i in range(self.k)]
                     for j in range(self.k)] for s in range(self.q)]

    def add(self, a, b):
        return tuple(x ^ y for x, y in zip(a, b))

    sub = add

    def neg(self, a):
        return a

    def scal(self, a, s):
        if s == 0:
            return self.ZERO
        if s == 1:
            return a
        m = self.mix[s]
        out = []
        for i in range(self.k):
            v = 0
            for j in range(self.k):
                if m[j][i]:
                    v ^= a[j]
            out.append(v)
        return tuple(out)

    def lincomb(self, cs, vs):
        # Bucket by coefficient value, XOR each bucket, then scale once -- so
        # the field multiplication happens q-1 times instead of once per term.
        if self.k == 2:
            # F_4 is the only prime-power field in the parameter sets, and this
            # is the hot loop of key generation; keeping the planes in locals
            # avoids building a tuple per term and runs several times faster
            # than the general path below.
            u0 = u1 = v0 = v1 = w0 = w1 = 0
            for i, c in enumerate(cs):
                if c:
                    p = vs[i]
                    if c == 1:
                        u0 ^= p[0]
                        u1 ^= p[1]
                    elif c == 2:
                        v0 ^= p[0]
                        v1 ^= p[1]
                    else:
                        w0 ^= p[0]
                        w1 ^= p[1]
            r = (u0, u1)
            if v0 or v1:
                r = self.add(r, self.scal((v0, v1), 2))
            if w0 or w1:
                r = self.add(r, self.scal((w0, w1), 3))
            return r
        buckets = [None] * self.q
        for i, c in enumerate(cs):
            if c:
                v = vs[i]
                acc = buckets[c]
                buckets[c] = v if acc is None else tuple(x ^ y for x, y in zip(acc, v))
        r = self.ZERO
        for s in range(1, self.q):
            if buckets[s] is not None:
                r = self.add(r, self.scal(buckets[s], s))
        return r

    def reduce(self, a):
        return a

    def polymul(self, a, b):
        k = self.k
        part = [[_clmul(a[i], b[j]) for j in range(k)] for i in range(k)]
        out = [0] * k
        for i in range(k):
            for j in range(k):
                p = part[i][j]
                if p:
                    # u^i * u^j lands on the planes named by mix[1<<i][j]
                    for t, bit in enumerate(self.mix[1 << i][j]):
                        if bit:
                            out[t] ^= p
        return tuple(out)

    def shift(self, a, n):
        return tuple(x << n for x in a)

    def rshift(self, a, n):
        return tuple(x >> n for x in a)

    def truncate(self, a, L):
        m = (1 << L) - 1
        return tuple(x & m for x in a)

    def coef(self, a, i):
        return sum(((a[j] >> i) & 1) << j for j in range(self.k))

    def to_list(self, a, L):
        return [self.coef(a, i) for i in range(L)]

    def from_list(self, cs):
        out = [0] * self.k
        for i, c in enumerate(cs):
            j = 0
            while c:
                if c & 1:
                    out[j] |= 1 << i
                c >>= 1
                j += 1
        return tuple(out)

    def is_zero(self, a):
        return not any(a)

    def degree(self, a):
        return max(x.bit_length() for x in a) - 1

    def support(self, a):
        d = self.degree(a)
        return [i for i in range(d + 1) if self.coef(a, i)]


class FpVec:
    """Vectors over F_p, p odd: coefficient i occupies lane i of B bits.

    Callers may let lanes run unreduced (values >= p) while accumulating; every
    routine that needs canonical coefficients calls `reduce` first.  B is sized
    so that a full polynomial product never overflows a lane.
    """

    def __init__(self, p, B):
        self.q = p
        self.B = B
        self.mask = (1 << B) - 1
        self.ZERO = 0

    def add(self, a, b):
        return a + b

    def neg(self, a):
        a = self.reduce(a)
        out, idx, B, p, m = 0, 0, self.B, self.q, self.mask
        while a:
            v = a & m
            if v:
                out |= (p - v) << (idx * B)
            a >>= B
            idx += 1
        return out

    def sub(self, a, b):
        return self.add(self.reduce(a), self.neg(b))

    def scal(self, a, s):
        return a * s

    def lincomb(self, cs, vs):
        r = 0
        for i, c in enumerate(cs):
            if c:
                r += c * vs[i]
        return self.reduce(r)

    def reduce(self, a):
        out, idx, B, p, m = 0, 0, self.B, self.q, self.mask
        while a:
            v = (a & m) % p
            if v:
                out |= v << (idx * B)
            a >>= B
            idx += 1
        return out

    def polymul(self, a, b):
        # One big-int multiply is the entire schoolbook product: lanes are wide
        # enough that nothing carries across.
        return self.reduce(a * b)

    def shift(self, a, k):
        return a << (k * self.B)

    def rshift(self, a, k):
        return a >> (k * self.B)

    def truncate(self, a, L):
        return a & ((1 << (L * self.B)) - 1)

    def coef(self, a, i):
        return (a >> (i * self.B)) & self.mask

    def to_list(self, a, L):
        B, m = self.B, self.mask
        return [(a >> (i * B)) & m for i in range(L)]

    def from_list(self, cs):
        r, B = 0, self.B
        for i, c in enumerate(cs):
            if c:
                r |= c << (i * B)
        return r

    def is_zero(self, a):
        return self.reduce(a) == 0

    def degree(self, a):
        a = self.reduce(a)
        return -1 if a == 0 else (a.bit_length() - 1) // self.B

    def support(self, a):
        return [i for i, c in enumerate(self.to_list(a, self.degree(a) + 1)) if c]


def make_vec(fq, n):
    """Pick the packed backend for F_q and a working length of about n."""
    q = fq.q
    if q == 2:
        return F2Vec()
    if fq.p == 2:
        return F2kVec(fq)
    # Lane width must hold n*(q-1)^2 (a full product) with room to spare.
    B = (n * (q - 1) ** 2).bit_length() + 2
    return FpVec(q, max(B, 8))


# ------------------------------------------------------- polynomials over Fq
#
# Only needed to decide irreducibility of the field polynomial; K-coefficient
# polynomials (the HFE central map) live in poly.py instead.


def _poly_rem(V, fq, a, b):
    """Remainder of a mod b, both packed vectors over F_q, b nonzero."""
    db = V.degree(b)
    binv = fq.inv(V.coef(V.reduce(b), db))
    a = V.reduce(a)
    while True:
        da = V.degree(a)
        if da < db:
            return a
        c = fq.mul(V.coef(a, da), binv)
        a = V.reduce(V.sub(a, V.shift(V.scal(b, c), da - db)))


def _poly_gcd(V, fq, a, b):
    a, b = V.reduce(a), V.reduce(b)
    while not V.is_zero(b):
        a, b = b, _poly_rem(V, fq, a, b)
    return a


# -------------------------------------------------------- the extension field


class ExtField:
    """K = F_q[t] / (f) with f monic irreducible of degree n.

    Elements are packed vectors of n coefficients over F_q; alpha denotes the
    class of t, so {1, alpha, ..., alpha^(n-1)} is the working basis and
    `to_list` gives exactly the F_q-coordinates the scheme's linear algebra
    operates on.
    """

    def __init__(self, q, n, fpoly=None, verify=True):
        self.q, self.n = q, n
        self.fq = Fq(q)
        self.V = make_vec(self.fq, 2 * n)
        self.fpoly = list(fpoly) if fpoly is not None else find_irreducible(q, n)
        if len(self.fpoly) != n:
            raise ValueError("field polynomial must have exactly n coefficients")
        self._build_reduction()
        if verify and not self.is_irreducible():
            raise ValueError("field polynomial is not irreducible")
        self._frob_cache = {}

    # -- setup ------------------------------------------------------------
    def _build_reduction(self):
        V, n = self.V, self.n
        # t^n = -(c_0 + c_1 t + ... + c_{n-1} t^{n-1})
        self.tn = V.reduce(V.neg(V.from_list(self.fpoly)))
        self.ZERO = V.ZERO
        self.ONE = V.from_list([1])
        self.ALPHA = V.from_list([0, 1]) if n > 1 else V.from_list([0])
        # RED[i] = t^(n+i) mod f, for i < n-1: everything reduction needs.
        red = [self.tn]
        for _ in range(n - 2):
            red.append(self.mul_t(red[-1]))
        self.RED = red

    def mul_t(self, a):
        """a * alpha, reduced."""
        V, n = self.V, self.n
        c = V.coef(a, n - 1)
        s = V.shift(V.truncate(a, n - 1), 1)
        if c:
            s = V.add(s, V.scal(self.tn, c))
        return V.reduce(s)

    # -- arithmetic -------------------------------------------------------
    def reduce_wide(self, w):
        """Reduce a product of degree < 2n-1 into K."""
        V, n = self.V, self.n
        w = V.reduce(w)
        hi = V.to_list(V.rshift(w, n), n - 1)
        lo = V.truncate(w, n)
        if any(hi):
            lo = V.add(lo, V.lincomb(hi, self.RED))
        return V.reduce(lo)

    def mul(self, a, b):
        return self.reduce_wide(self.V.polymul(a, b))

    def add(self, a, b):
        return self.V.reduce(self.V.add(a, b))

    def sub(self, a, b):
        return self.V.reduce(self.V.sub(a, b))

    def neg(self, a):
        return self.V.neg(a)

    def scal(self, a, s):
        return self.V.reduce(self.V.scal(a, s))

    def is_zero(self, a):
        return self.V.is_zero(a)

    def eq(self, a, b):
        return self.V.reduce(a) == self.V.reduce(b)

    def pow(self, a, e):
        r, b = self.ONE, a
        while e:
            if e & 1:
                r = self.mul(r, b)
            e >>= 1
            if e:
                b = self.mul(b, b)
        return r

    def inv(self, a):
        if self.is_zero(a):
            raise ZeroDivisionError("0 has no inverse in K")
        return self.pow(a, self.q ** self.n - 2)

    def coords(self, a):
        return self.V.to_list(self.V.reduce(a), self.n)

    def from_coords(self, cs):
        return self.V.reduce(self.V.from_list(list(cs)))

    def random(self, xof):
        from .symmetric import sample_fq
        return self.from_coords(sample_fq(xof, self.n, self.q))

    def random_nonzero(self, xof):
        while True:
            a = self.random(xof)
            if not self.is_zero(a):
                return a

    # -- Frobenius --------------------------------------------------------
    def frob_table(self, k):
        """[ (alpha^i)^(q^k) for i < n ].

        Since every F_q-coordinate is fixed by x -> x^q, raising an arbitrary
        element to the q^k is just a linear combination against this table.
        """
        if k in self._frob_cache:
            return self._frob_cache[k]
        beta = self.ALPHA
        for _ in range(k % self.n):
            beta = self.pow(beta, self.q)
        tbl, cur = [], self.ONE
        for _ in range(self.n):
            tbl.append(cur)
            cur = self.mul(cur, beta)
        self._frob_cache[k] = tbl
        return tbl

    def frob(self, a, k):
        """a^(q^k)."""
        k %= self.n
        if k == 0:
            return self.V.reduce(a)
        return self.V.lincomb(self.coords(a), self.frob_table(k))

    # -- irreducibility ---------------------------------------------------
    def is_irreducible(self):
        """Rabin's test for f, using the q-power Frobenius on F_q[t]/(f).

        On the quotient ring, h -> h^q is F_q-linear (the coefficients lie in
        F_q and so are fixed), hence one application is a single linear
        combination against the precomputed powers t^(q i) mod f.  That makes
        each of the n iterations O(n) rather than a modular exponentiation.
        """
        V, n, q = self.V, self.n, self.q
        if n == 1:
            return True
        # TQ[i] = t^(q i) mod f
        step = self.ONE
        for _ in range(q):
            step = self.mul_t(step)          # step = t^q mod f
        TQ, cur = [], self.ONE
        for _ in range(n):
            TQ.append(cur)
            cur = self.mul(cur, step)

        def phi(h):
            return V.lincomb(V.to_list(h, n), TQ)

        primes, m = set(), n
        d = 2
        while d * d <= m:
            if m % d == 0:
                primes.add(d)
                while m % d == 0:
                    m //= d
            d += 1
        if m > 1:
            primes.add(m)
        checkpoints = {n // p for p in primes}

        cur = self.ALPHA
        for k in range(1, n + 1):
            cur = phi(cur)
            if k in checkpoints:
                g = _poly_gcd(V, self.fq, V.sub(cur, self.ALPHA),
                              V.add(V.from_list(self.fpoly), V.shift(self.ONE, n)))
                if V.degree(g) != 0:
                    return False
        return self.V.reduce(cur) == self.V.reduce(self.ALPHA)


def find_irreducible(q, n):
    """The canonical field polynomial for (q, n).

    Deterministic: candidates f = t^n + sum c_i t^i are enumerated by the
    integer whose base-q digits are (c_0, ..., c_{n-1}), and the first
    irreducible one wins.  Candidates with c_0 = 0 are skipped, being
    divisible by t.
    """
    key = 0
    while True:
        key += 1
        cs, v = [], key
        for _ in range(n):
            cs.append(v % q)
            v //= q
        if v:                                  # ran out of room: no such f
            raise AssertionError("no irreducible found for q=%d n=%d" % (q, n))
        if cs[0] == 0:
            continue
        cand = ExtField(q, n, fpoly=cs, verify=False)
        if cand.is_irreducible():
            return cs


# ------------------------------------------------- linear-combination helper
#
# The public-key loop evaluates the same shape over and over: take the F_q
# coordinates of one vector and combine a table of vectors against them.  Over
# F_2 only the *positions* of the nonzero coordinates matter, which halves the
# work, so `prep` hands back whatever form the matching `lincomb_prep` wants.
# The prepared form is backend-specific but space-independent, so coordinates
# taken in K can drive a combination in F_q^m.


def prep(V, a, L):
    if isinstance(V, F2Vec):
        return V.support(V.truncate(a, L))
    return V.to_list(V.reduce(a), L)


def lincomb_prep(W, prepped, tbl):
    if isinstance(W, F2Vec):
        r = 0
        for c in prepped:
            r ^= tbl[c]
        return r
    return W.lincomb(prepped, tbl)

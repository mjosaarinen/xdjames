"""Test suite.  Run with:  python3 -m unittest discover -s tests -v"""

import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from djames import codec, kat, linalg, params, scheme  # noqa: E402
from djames import poly as P  # noqa: E402
from djames.ff import ExtField, Fq, make_vec  # noqa: E402
import math  # noqa: E402
from djames.symmetric import XOF, sample_fq  # noqa: E402

FIELDS = (2, 4, 5, 13, 23)
KATDIR = os.path.join(os.path.dirname(__file__), "..", "kat")


class TestParameters(unittest.TestCase):
    def test_confirmed_parameters(self):
        for name in params.names():
            p = params.get(name)
            if p.q == 4:
                self.assertEqual(p.monomials, [(0, 1), (0, 2)], name)
                self.assertEqual(p.tag().decode().rsplit("/", 1)[-1],
                                 "mon0-1.0-2", name)

        p = params.get("d-james-256-q4")
        self.assertEqual((p.m, p.a, p.n), (170, 53, 223))

    def test_monomials_are_distinct(self):
        with self.assertRaises(AssertionError):
            params.Params("duplicate", "james", 5, 8, 2, 1, 6,
                          [(0, 0), (0, 0)])


class TestBaseField(unittest.TestCase):
    def test_axioms(self):
        for q in FIELDS:
            f = Fq(q)
            for a in range(q):
                self.assertEqual(f.add(a, f.neg(a)), 0)
                if a:
                    self.assertEqual(f.mul(a, f.inv(a)), 1)
                for b in range(q):
                    for c in range(q):
                        self.assertEqual(f.mul(a, f.add(b, c)),
                                         f.add(f.mul(a, b), f.mul(a, c)))


class TestPackedVectors(unittest.TestCase):
    def test_polymul_matches_schoolbook(self):
        rng = random.Random(1)
        for q in FIELDS:
            fq, V = Fq(q), make_vec(Fq(q), 60)
            for _ in range(25):
                L = rng.randint(1, 30)
                a = [rng.randrange(q) for _ in range(L)]
                b = [rng.randrange(q) for _ in range(L)]
                ref = [0] * (2 * L - 1)
                for i, x in enumerate(a):
                    for j, y in enumerate(b):
                        ref[i + j] = fq.add(ref[i + j], fq.mul(x, y))
                got = V.to_list(V.reduce(V.polymul(V.from_list(a),
                                                   V.from_list(b))), 2 * L - 1)
                self.assertEqual(got, ref, (q, a, b))


class TestExtField(unittest.TestCase):
    def test_field_and_frobenius(self):
        rng = random.Random(2)
        for q, n in ((2, 16), (2, 21), (4, 7), (5, 9), (13, 6), (23, 5)):
            K = ExtField(q, n)
            rnd = lambda: K.from_coords([rng.randrange(q) for _ in range(n)])
            for _ in range(20):
                a, b, c = rnd(), rnd(), rnd()
                self.assertTrue(K.eq(K.mul(a, K.add(b, c)),
                                     K.add(K.mul(a, b), K.mul(a, c))))
                if not K.is_zero(a):
                    self.assertTrue(K.eq(K.mul(a, K.inv(a)), K.ONE))
                self.assertTrue(K.eq(K.frob(a, n), a))
                k = rng.randrange(1, n)
                ref = a
                for _ in range(k):
                    ref = K.pow(ref, q)
                self.assertTrue(K.eq(K.frob(a, k), ref))

    def test_field_polynomials_are_irreducible(self):
        for name in params.names():
            p = params.get(name)
            K = ExtField(p.q, p.n, fpoly=p.fpoly, verify=False)
            self.assertTrue(K.is_irreducible(), name)


class TestRootFinding(unittest.TestCase):
    def test_planted_roots_are_found(self):
        rng = random.Random(3)
        x = XOF(b"unit-test")
        for q, n, D in ((2, 17, 5), (4, 9, 17), (5, 11, 6), (13, 7, 14), (23, 6, 24)):
            K = ExtField(q, n)
            rnd = lambda: K.from_coords([rng.randrange(q) for _ in range(n)])
            for _ in range(8):
                planted = [rnd() for _ in range(rng.randint(0, 3))]
                f = [K.ONE]
                for r in planted:
                    f = P.mul(K, f, [K.neg(r), K.ONE])
                while P.deg(f) < D:
                    f = P.mul(K, f, [rnd(), rnd(), K.ONE])
                found = P.roots(K, f, x)
                for r in found:
                    self.assertTrue(K.is_zero(P.evaluate(K, f, r)))
                got = {K.V.reduce(r) for r in found}
                for r in planted:
                    self.assertIn(K.V.reduce(r), got)


class TestLinalg(unittest.TestCase):
    def test_inverse(self):
        for q in FIELDS:
            fq, n = Fq(q), 12
            M, Mi = linalg.random_invertible(fq, n, XOF(b"la", bytes([q])))
            for i in range(n):
                e = [1 if k == i else 0 for k in range(n)]
                self.assertEqual(linalg.vec_mat(fq, linalg.vec_mat(fq, e, M), Mi), e)
        for q in FIELDS:
            fq = Fq(q)
            Z = linalg.random_full_rank(fq, 3, 12, XOF(b"z", bytes([q])))
            self.assertEqual(linalg.rank(fq, Z, 12), 3)


class TestCodec(unittest.TestCase):
    def test_grouped_round_trip(self):
        rng = random.Random(4)
        for q in FIELDS:
            for _ in range(15):
                n = rng.randint(1, 300)
                d = [rng.randrange(q) for _ in range(n)]
                b = codec.digits_to_bytes(d, q)
                self.assertEqual(len(b), codec.packed_len(n, q))
                self.assertEqual(codec.bytes_to_digits(b, n, q), d)

    def test_vector_round_trip(self):
        rng = random.Random(8)
        for q in FIELDS:
            for _ in range(20):
                n = rng.randint(1, 400)
                d = [rng.randrange(q) for _ in range(n)]
                b = codec.encode_vector(d, q)
                self.assertEqual(len(b), codec.vec_bytes(n, q))
                self.assertEqual(codec.decode_vector(b, n, q), d)

    def test_vector_length_is_information_theoretic(self):
        """The signature encoding must not waste a byte."""
        for q in FIELDS:
            for n in (1, 7, 74, 94, 171, 189, 390, 578):
                self.assertEqual(codec.vec_bytes(n, q),
                                 math.ceil(n * math.log2(q) / 8))

    def test_vector_encoding_is_canonical(self):
        """Every byte string that is not the unique encoding must be rejected.

        Otherwise signatures are malleable: over F_2 the unused high bits of
        the last byte are free, and over odd F_q one may add q^n.
        """
        rng = random.Random(5)
        for q in FIELDS:
            for _ in range(10):
                n = rng.randint(4, 120)
                d = [rng.randrange(q) for _ in range(n)]
                b = codec.encode_vector(d, q)
                v = int.from_bytes(b, "little")
                for extra in (q ** n, 2 * q ** n):
                    if v + extra < 256 ** len(b):
                        alt = (v + extra).to_bytes(len(b), "little")
                        self.assertRaises(ValueError,
                                          codec.decode_vector, alt, n, q)
                self.assertRaises(ValueError, codec.decode_vector, b + b"\x00", n, q)
                self.assertRaises(ValueError, codec.decode_vector, b[:-1], n, q)

    def test_grouped_encoding_is_canonical(self):
        rng = random.Random(6)
        for q in FIELDS:
            g, B = codec.pack_shape(q)
            n = 3 * g + 3
            d = [rng.randrange(q) for _ in range(n)]
            b = bytearray(codec.digits_to_bytes(d, q))
            self.assertRaises(ValueError, codec.bytes_to_digits,
                              bytes(b) + b"\x00", n, q)
            if q & (q - 1):                   # power of two: every byte is valid
                v = int.from_bytes(bytes(b[:B]), "little") + q ** g
                if v < 256 ** B:
                    b[:B] = v.to_bytes(B, "little")
                    self.assertRaises(ValueError, codec.bytes_to_digits,
                                      bytes(b), n, q)

class TestScheme(unittest.TestCase):
    """Toy parameters only -- correctness, not security."""

    def test_seed_length_is_enforced(self):
        p = params.get("toy-d-james-q2")
        for bad in (b"", b"short", b"x" * 31, b"x" * 33):
            self.assertRaises(ValueError, scheme.keygen, p, bad)
        _, sk = scheme.keygen(p, b"x" * scheme.SEED_BYTES)
        self.assertEqual(len(sk.to_bytes()), scheme.SEED_BYTES)

    def _one(self, name):
        p = params.get(name)
        pk, sk = scheme.keygen(p, b"unit-test-seed--" * 2)
        raw = pk.to_bytes()

        # deterministic key generation
        pk2, _ = scheme.keygen(p, b"unit-test-seed--" * 2)
        self.assertEqual(pk2.to_bytes(), raw)
        pk3, _ = scheme.keygen(p, b"different-seed--" * 2)
        self.assertNotEqual(pk3.to_bytes(), raw)

        # public key survives serialization
        self.assertEqual(scheme.PublicKey.from_bytes(p, raw).to_bytes(), raw)
        loaded = scheme.PublicKey.from_bytes(p, raw)

        for i in range(2):
            msg = b"message number %d" % i
            sig = scheme.sign(p, sk, msg)
            self.assertEqual(len(sig), p.sig_bytes)
            self.assertEqual(len(sig), math.ceil(p.sig_bits / 8))
            self.assertEqual(scheme.sign(p, sk, msg), sig)      # deterministic
            self.assertTrue(scheme.verify(p, pk, msg, sig))
            self.assertTrue(scheme.verify(p, loaded, msg, sig))
            self.assertFalse(scheme.verify(p, pk, msg + b"!", sig))
            bad = bytearray(sig)
            bad[0] ^= 1
            self.assertFalse(scheme.verify(p, pk, msg, bytes(bad)))
            self.assertFalse(scheme.verify(p, pk, msg, b"\x00" * len(sig)))
            self.assertFalse(scheme.verify(p, pk, msg, sig + b"\x00"))
            self.assertFalse(scheme.verify(p, pk, msg, sig[:-1]))

            # Non-malleability: the signature has exactly one valid encoding.
            v = int.from_bytes(sig, "little")
            for extra in (p.q ** p.n, 2 * p.q ** p.n):
                if v + extra < 256 ** len(sig):
                    alt = (v + extra).to_bytes(len(sig), "little")
                    self.assertFalse(scheme.verify(p, pk, msg, alt))
            for bit in range(8):
                alt = bytearray(sig)
                alt[-1] ^= 1 << bit
                if bytes(alt) != sig:
                    self.assertFalse(scheme.verify(p, pk, msg, bytes(alt)),
                                     "%s: alternate encoding accepted" % name)


def _add_scheme_tests():
    for name in params.names("toy-"):
        setattr(TestScheme, "test_" + name.replace("-", "_"),
                lambda self, n=name: self._one(n))


_add_scheme_tests()


class TestKAT(unittest.TestCase):
    def test_vectors(self):
        found = False
        for fn in sorted(os.listdir(KATDIR)) if os.path.isdir(KATDIR) else []:
            if not fn.endswith(".json"):
                continue
            # The real-parameter vectors are long-running; opt in explicitly.
            if fn != "toy.json" and not os.environ.get("DJAMES_FULL_KAT"):
                continue
            for vec in kat.load(os.path.join(KATDIR, fn)):
                found = True
                self.assertEqual(kat.check(vec), [], "%s / %s" % (fn, vec["name"]))
        if not found:
            self.skipTest("no KAT files (run tools/gen_kat.py toy)")


if __name__ == "__main__":
    unittest.main()


class TestPublicKeyMatchesCentralMap(unittest.TestCase):
    """The strongest available check on key generation.

    Signing never touches the public key and verification never touches the
    secret key, so a sign/verify round-trip already says a lot -- but it only
    ever checks that the public system evaluates to *zero*.  Here we evaluate
    the secret central map directly at uniformly random (a, b), project it
    through T, and require it to agree with the public key coefficient-for-
    coefficient.  That pins every published coefficient, not just the ones a
    signature happens to exercise.
    """

    @staticmethod
    def _central(p, sk, a, b):
        K, fq = sk.K, sk.K.fq
        x = linalg.vec_mat(fq, a, sk.S)                  # x = a S
        X = K.from_coords(x)
        z = linalg.mat_vec(fq, sk.MZ, x)
        acc = K.ZERO
        for (i, j), lm in zip(p.monomials, sk.lam):      # HFE core
            acc = K.add(acc, K.mul(lm, K.mul(K.frob(X, i), K.frob(X, j))))
        if p.ny is not None:                             # Dragon
            for k in range(p.d + 1):
                Lk = K.ZERO
                for t in range(p.ny):
                    if b[t]:
                        Lk = K.add(Lk, K.scal(sk.Lam[k][t], b[t]))
                acc = K.add(acc, K.mul(Lk, K.frob(X, k)))
        for k in range(p.d):                             # IP, bilinear
            s = K.ZERO
            for t in range(p.r):
                if z[t]:
                    s = K.add(s, K.scal(sk.Mbil[k][t], z[t]))
            acc = K.add(acc, K.mul(s, K.frob(X, k)))
        for i in range(p.r):                             # IP, quadratic
            for j in range(i, p.r):
                c = fq.mul(z[i], z[j])
                if c:
                    acc = K.add(acc, K.scal(sk.G[i][j], c))
        co = K.coords(acc)
        out = []
        for k in range(p.m):                             # project, then minus
            v = 0
            for c in range(p.n):
                if sk.T[k][c] and co[c]:
                    v = fq.add(v, fq.mul(sk.T[k][c], co[c]))
            out.append(v)
        return out

    def _one(self, name):
        p = params.get(name)
        pk, sk = scheme.keygen(p, b"central-map-check" + b"\x00" * 15)
        x = XOF(b"central-map-inputs", name.encode())
        for _ in range(4):
            a = sample_fq(x, p.n, p.q)
            b = sample_fq(x, p.ny, p.q) if p.ny else None
            want = self._central(p, sk, a, b)
            got = pk.W.to_list(scheme._evaluate(p, pk, a, b), p.m)
            self.assertEqual(got, want, name)


for _n in params.names("toy-"):
    setattr(TestPublicKeyMatchesCentralMap, "test_" + _n.replace("-", "_"),
            lambda self, n=_n: self._one(n))

"""Known-answer tests: generation and checking.

A vector pins the public key by its SHA3-256 digest (the keys run to megabytes)
plus its exact length, and pins signatures in full.  Because keygen and signing
are deterministic functions of the seed and the message, a vector either
reproduces bit for bit or the implementation changed.
"""

import hashlib
import json

from . import params as _params
from .scheme import keygen, sign, verify, PublicKey


def _digest(b):
    return hashlib.sha3_256(b).hexdigest()


MESSAGES = [b"", b"abc", bytes(range(32)),
            b"D-James known-answer test vector"]


def generate(name, seed=None, nmsg=len(MESSAGES), sign_vectors=True):
    P = _params.get(name)
    if seed is None:
        seed = hashlib.shake_256(name.encode()).digest(32)
    pk, sk = keygen(P, seed)
    raw = pk.to_bytes()
    out = {
        "name": name,
        "scheme": P.scheme,
        "params": {"q": P.q, "n": P.n, "m": P.m, "a": P.a, "ny": P.ny,
                   "r": P.r, "D": P.D, "d": P.d, "monomials": P.monomials},
        "field_poly_support": [i for i, c in enumerate(P.fpoly) if c],
        "seed": seed.hex(),
        "pk_len": len(raw),
        "pk_sha3_256": _digest(raw),
        "sig_bits": P.sig_bits,
        "signatures": [],
    }
    if sign_vectors:
        for msg in MESSAGES[:nmsg]:
            sig = sign(P, sk, msg)
            assert verify(P, pk, msg, sig)
            out["signatures"].append({"msg": msg.hex(), "sig": sig.hex()})
    return out


def check(vec, resign=True):
    """Re-derive a vector and compare.  Returns a list of mismatch strings."""
    P = _params.get(vec["name"])
    bad = []
    seed = bytes.fromhex(vec["seed"])
    pk, sk = keygen(P, seed)
    raw = pk.to_bytes()
    if len(raw) != vec["pk_len"]:
        bad.append("pk_len %d != %d" % (len(raw), vec["pk_len"]))
    if _digest(raw) != vec["pk_sha3_256"]:
        bad.append("pk digest mismatch")
    for i, s in enumerate(vec["signatures"]):
        msg, want = bytes.fromhex(s["msg"]), bytes.fromhex(s["sig"])
        # The recorded signature must verify even if we do not re-sign.
        if not verify(P, pk, msg, want):
            bad.append("recorded signature %d does not verify" % i)
        if resign:
            got = sign(P, sk, msg)
            if got != want:
                bad.append("signature %d: %s != %s" % (i, got.hex(), want.hex()))
    return bad


def load(path):
    with open(path) as f:
        return json.load(f)


def save(path, vectors):
    with open(path, "w") as f:
        json.dump(vectors, f, indent=1, sort_keys=True)
        f.write("\n")

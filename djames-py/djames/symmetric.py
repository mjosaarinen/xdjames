"""SHAKE-based symmetric layer.

Every byte of randomness in this implementation -- secret key material, the
message digest, the auxiliary randomness inside equal-degree splitting -- is
squeezed from SHAKE256.  There is no other entropy source, so `keygen(seed)`
and `sign(sk, msg)` are bit-deterministic functions of their inputs, which is
what makes reproducible test vectors possible.

Domain separation: every stream is opened with a short ASCII label so that
two different uses of the same seed never produce the same stream.
"""

import hashlib

# Every label used anywhere in the scheme, kept together so the domains are
# easy to audit for collisions.
DOM_KEY = b"D-James/v1/keygen"        # expand the master seed
DOM_PRF = b"D-James/v1/prf"           # derive the per-key signing PRF key
DOM_MSG = b"D-James/v1/msg"           # hash a message to F_q^{n_y}
DOM_PAD = b"D-James/v1/pad"           # James: the random minus-padding
DOM_EDF = b"D-James/v1/edf"           # equal-degree splitting randomness


class XOF:
    """A readable SHAKE256 stream.

    `hashlib` exposes SHAKE as a fixed-length digest rather than as a
    squeezable sponge, but `shake_256(x).digest(k)` is by definition the first
    k bytes of the same stream, so a longer digest always extends a shorter
    one.  Buffer doubling keeps the total re-hashing work linear.
    """

    __slots__ = ("_h", "_buf", "_pos")

    def __init__(self, *parts: bytes):
        h = hashlib.shake_256()
        for p in parts:
            # Length-prefix each part so that (b"ab", b"c") and (b"a", b"bc")
            # are distinct inputs.
            h.update(len(p).to_bytes(4, "little"))
            h.update(p)
        self._h = h
        self._buf = b""
        self._pos = 0

    def read(self, nbytes: int) -> bytes:
        end = self._pos + nbytes
        if end > len(self._buf):
            self._buf = self._h.digest(max(64, 2 * len(self._buf), end))
        out = self._buf[self._pos:end]
        self._pos = end
        return out

    def read_int(self, nbytes: int) -> int:
        return int.from_bytes(self.read(nbytes), "little")


def _bits_per(q: int) -> int:
    """log2(q) when q is a power of two, else 0."""
    return q.bit_length() - 1 if q & (q - 1) == 0 else 0


def sample_fq(xof: XOF, count: int, q: int) -> list:
    """`count` uniform elements of F_q, as ints in [0, q).

    For q a power of two the field elements are read straight off the bit
    stream, consuming exactly ceil(count * log2 q / 8) bytes.  Otherwise we
    rejection-sample one byte at a time from the largest multiple of q below
    256, which keeps the output exactly uniform and (for q <= 23) wastes
    under 10% of the stream.

    Both paths consume a precisely defined number of stream bytes, which
    matters because callers keep drawing from the same XOF afterwards.
    """
    k = _bits_per(q)
    if k:
        out = []
        need = (count * k + 7) // 8
        acc = int.from_bytes(xof.read(need), "little")
        mask = q - 1
        for _ in range(count):
            out.append(acc & mask)
            acc >>= k
        return out

    # One byte at a time.  Reading in larger chunks would give the same
    # digits but would advance the stream by a chunk-size-dependent amount,
    # so every later draw from the same XOF would diverge between two
    # implementations that chunk differently.  Byte-at-a-time is the only
    # choice that needs no further specification.
    limit = (256 // q) * q
    out = []
    while len(out) < count:
        b = xof.read(1)[0]
        if b < limit:
            out.append(b % q)
    return out


def hash_message(msg: bytes, salt: int, count: int, q: int) -> list:
    """Hash (msg, salt) to `count` elements of F_q.

    `salt` is the D-James signing counter.  It is *not* transmitted: the
    verifier re-derives it by trying salt = 0, 1, 2, ... (see the paper's
    footnote 5, which trades verification time for a shorter signature).
    """
    x = XOF(DOM_MSG, msg, salt.to_bytes(8, "little"), q.to_bytes(4, "little"),
            count.to_bytes(4, "little"))
    return sample_fq(x, count, q)

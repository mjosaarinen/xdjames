"""Deterministic byte encodings.

Two encoders, because signatures and public keys want different things.

`encode_vector` / `decode_vector` -- used for **signatures**.  The whole vector
is one base-q integer, so the encoding is exactly ceil(n log2 q / 8) bytes and,
crucially, *canonical*: `decode_vector` rejects any byte string that is not the
unique encoding of its digit vector.  Without that check a signature would be
malleable -- over F_2 the unused high bits of the final byte are free, and over
odd F_q one can add q^g to a group -- so distinct byte strings would verify
against the same message.

`DigitPacker` / `bytes_to_digits` -- used for **public keys**, which run to
millions of digits and would make a single base conversion quadratic.  Digits
are packed in groups of g into B bytes, with (g, B) chosen by the search below
to land within about 1% of the information-theoretic size; the trailing partial
group is packed exactly, not padded out to a full one.  Decoding validates
every group, so this encoding is canonical too.
"""

from math import gcd


def pack_shape(q):
    """(digits per group, bytes per group)."""
    if q & (q - 1) == 0:
        k = q.bit_length() - 1
        l = (k * 8) // gcd(k, 8)
        return l // k, l // 8
    best = None
    for B in range(1, 17):
        cap = 1 << (8 * B)
        g, v = 0, 1
        while v * q <= cap:
            v *= q
            g += 1
        if g and (best is None or g / B > best[0] / best[1]):
            best = (g, B)
    return best


def vec_bytes(count, q):
    """Exact byte length of `count` F_q digits: ceil(count * log2 q / 8)."""
    if count <= 0:
        return 0
    return ((q ** count - 1).bit_length() + 7) // 8


def encode_vector(digits, q):
    """The canonical encoding of a digit vector: one base-q integer."""
    v = 0
    for d in reversed(digits):
        v = v * q + d
    return v.to_bytes(vec_bytes(len(digits), q), "little")


def decode_vector(data, count, q):
    """Inverse of encode_vector, rejecting every non-canonical encoding."""
    if len(data) != vec_bytes(count, q):
        raise ValueError("expected %d bytes, got %d"
                         % (vec_bytes(count, q), len(data)))
    v = int.from_bytes(data, "little")
    if v >= q ** count:
        raise ValueError("non-canonical encoding")
    out = []
    for _ in range(count):
        out.append(v % q)
        v //= q
    return out


def packed_len(count, q):
    """Length of the grouped encoding, with the final group packed exactly."""
    g, B = pack_shape(q)
    full, rem = divmod(count, g)
    return full * B + vec_bytes(rem, q)


class DigitPacker:
    """Append F_q digits, get bytes.

    Vectors arrive already packed (a Python int over F_2 / F_{2^k}, or lanes
    over F_p), so the power-of-two path never materialises a digit list: it
    shifts the packed value straight into a bit accumulator.
    """

    def __init__(self, q):
        self.q = q
        self.g, self.B = pack_shape(q)
        self.out = bytearray()
        self.pow2 = (q & (q - 1)) == 0
        self.k = q.bit_length() - 1 if self.pow2 else 0
        self.acc = 0
        self.nbits = 0
        self.pending = []

    def push_digits(self, digits):
        if self.pow2:
            for d in digits:
                self.acc |= d << self.nbits
                self.nbits += self.k
            self._flush_bits()
        else:
            self.pending.extend(digits)
            self._flush_groups()

    def push_packed(self, V, vec, L):
        """Append the L low coefficients of a packed vector."""
        if self.pow2 and self.k == 1:
            self.acc |= (V.truncate(vec, L)) << self.nbits
            self.nbits += L
            self._flush_bits()
        else:
            self.push_digits(V.to_list(V.reduce(vec), L))

    def _flush_bits(self):
        while self.nbits >= 64:
            self.out += (self.acc & 0xFFFFFFFFFFFFFFFF).to_bytes(8, "little")
            self.acc >>= 64
            self.nbits -= 64

    def _flush_groups(self):
        g, B, q = self.g, self.B, self.q
        p = self.pending
        i = 0
        while i + g <= len(p):
            v = 0
            for j in range(g - 1, -1, -1):
                v = v * q + p[i + j]
            self.out += v.to_bytes(B, "little")
            i += g
        del p[:i]

    def bytes(self):
        out = bytearray(self.out)
        if self.pow2:
            if self.nbits:
                out += self.acc.to_bytes((self.nbits + 7) // 8, "little")
        elif self.pending:
            # Pack the tail into exactly the bytes it needs.  Padding it out to
            # a full group is what used to cost up to 11 bytes per vector.
            out += encode_vector(self.pending, self.q)
        return bytes(out)


def digits_to_bytes(digits, q):
    p = DigitPacker(q)
    p.push_digits(list(digits))
    return p.bytes()


def bytes_to_digits(data, count, q):
    """Inverse of digits_to_bytes, rejecting non-canonical encodings."""
    if len(data) != packed_len(count, q):
        raise ValueError("expected %d bytes, got %d"
                         % (packed_len(count, q), len(data)))
    g, B = pack_shape(q)
    out, pos = [], 0
    while len(out) < count:
        take = min(g, count - len(out))
        nb = B if take == g else vec_bytes(take, q)
        v = int.from_bytes(data[pos:pos + nb], "little")
        pos += nb
        if v >= q ** take:
            raise ValueError("non-canonical group at digit %d" % len(out))
        for _ in range(take):
            out.append(v % q)
            v //= q
    return out

//! The symmetric layer: domain-separated SHAKE256 streams and F_q sampling.
//!
//! Every byte of randomness in the scheme comes from here, so `keygen` and
//! `sign` are deterministic functions of their inputs. See `d-james-spec.md`
//! §4; the byte counts below are normative, because callers keep drawing from
//! the same stream afterwards.

use crate::keccak::Shake256;
use alloc::vec::Vec;

pub const DOM_KEY: &[u8] = b"D-James/v1/keygen";
pub const DOM_PRF: &[u8] = b"D-James/v1/prf";
pub const DOM_MSG: &[u8] = b"D-James/v1/msg";
pub const DOM_PAD: &[u8] = b"D-James/v1/pad";
pub const DOM_EDF: &[u8] = b"D-James/v1/edf";

/// A readable SHAKE256 stream over length-prefixed parts.
///
/// Each part is prefixed with its length as four little-endian bytes, so
/// `("ab", "c")` and `("a", "bc")` are distinct inputs.
pub struct Xof {
    s: Shake256,
    buf: [u8; 136],
    pos: usize,
    len: usize,
}

impl Xof {
    pub fn new(parts: &[&[u8]]) -> Self {
        let mut s = Shake256::new();
        for p in parts {
            s.absorb(&(p.len() as u32).to_le_bytes());
            s.absorb(p);
        }
        Xof {
            s,
            buf: [0u8; 136],
            pos: 0,
            len: 0,
        }
    }

    pub fn read(&mut self, out: &mut [u8]) {
        let mut off = 0;
        while off < out.len() {
            if self.pos == self.len {
                self.s.squeeze(&mut self.buf);
                self.pos = 0;
                self.len = self.buf.len();
            }
            let n = core::cmp::min(out.len() - off, self.len - self.pos);
            out[off..off + n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
            self.pos += n;
            off += n;
        }
    }

    #[inline]
    pub fn byte(&mut self) -> u8 {
        if self.pos == self.len {
            self.s.squeeze(&mut self.buf);
            self.pos = 0;
            self.len = self.buf.len();
        }
        let b = self.buf[self.pos];
        self.pos += 1;
        b
    }
}

/// `count` uniform elements of F_q, as bytes in `[0, q)`.
///
/// For `q` a power of two the digits are read straight off the bit stream,
/// consuming exactly `ceil(count * log2 q / 8)` bytes. Otherwise bytes are
/// rejection-sampled **one at a time** against the largest multiple of `q`
/// below 256. Reading the rejection path in larger chunks would give the same
/// digits but leave the stream at a different offset, so byte-at-a-time is
/// normative rather than incidental.
pub fn sample_fq(xof: &mut Xof, count: usize, q: u32) -> Vec<u8> {
    let mut out = Vec::with_capacity(count);
    if q.is_power_of_two() {
        let k = q.trailing_zeros() as usize;
        let need = (count * k + 7) / 8;
        let mut buf = alloc::vec![0u8; need];
        xof.read(&mut buf);
        let mask = (q - 1) as u8;
        let mut bit = 0usize;
        for _ in 0..count {
            let byte = bit / 8;
            let sh = bit % 8;
            // k is 1 or 2 and 8 is a multiple of both, so a digit never
            // straddles a byte boundary.
            out.push((buf[byte] >> sh) & mask);
            bit += k;
        }
        crate::ct::zeroize8(&mut buf);
    } else {
        let limit = (256 / q) * q;
        while out.len() < count {
            let b = xof.byte() as u32;
            if b < limit {
                out.push((b % q) as u8);
            }
        }
    }
    out
}

/// Hash `(msg, salt)` to `count` elements of F_q.
pub fn hash_to_fq(msg: &[u8], salt: u64, count: usize, q: u32) -> Vec<u8> {
    let mut x = Xof::new(&[
        DOM_MSG,
        msg,
        &salt.to_le_bytes(),
        &q.to_le_bytes(),
        &(count as u32).to_le_bytes(),
    ]);
    sample_fq(&mut x, count, q)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cross-check against the Python reference implementation.  These values
    /// were produced by `djames-py`; if either side's XOF framing or sampler
    /// drifts, every derived key would silently diverge, so they are pinned.
    #[test]
    fn matches_python_reference() {
        let mut x = Xof::new(&[DOM_KEY, b"seed", b"tag"]);
        let mut b = [0u8; 16];
        x.read(&mut b);
        assert_eq!(
            b,
            [
                0xb5, 0x75, 0xbc, 0x69, 0x61, 0x47, 0xc5, 0x81, 0x9c, 0x94, 0x92, 0x7e, 0x5f, 0x35,
                0x89, 0x6c
            ]
        );

        let cases: [(u32, [u8; 12]); 5] = [
            (2, [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1]),
            (4, [3, 0, 2, 1, 1, 2, 3, 0, 1, 3, 0, 1]),
            (5, [4, 2, 2, 2, 2, 1, 2, 1, 0, 4, 3, 2]),
            (13, [8, 5, 12, 4, 10, 6, 6, 4, 0, 12, 8, 2]),
            (23, [7, 11, 8, 17, 17, 5, 9, 6, 3, 16, 7, 1]),
        ];
        for (q, want) in cases {
            let mut x = Xof::new(&[b"a"]);
            assert_eq!(sample_fq(&mut x, 12, q), want.to_vec(), "q={q}");
        }

        assert_eq!(
            hash_to_fq(b"hi", 0, 8, 5),
            alloc::vec![4, 2, 1, 0, 1, 0, 1, 3]
        );
        assert_eq!(
            hash_to_fq(b"hi", 3, 8, 23),
            alloc::vec![16, 20, 8, 13, 8, 18, 20, 3]
        );
    }

    #[test]
    fn read_chunking_is_irrelevant() {
        let mut a = [0u8; 300];
        Xof::new(&[b"x", b"y"]).read(&mut a);
        let mut x = Xof::new(&[b"x", b"y"]);
        let mut b = [0u8; 300];
        let mut off = 0;
        for n in [1usize, 5, 130, 136, 27, 101] {
            let n = n.min(300 - off);
            x.read(&mut b[off..off + n]);
            off += n;
        }
        assert_eq!(off, 300);
        assert_eq!(a, b);
    }
}

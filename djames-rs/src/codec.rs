//! Deterministic, canonical byte encodings (`d-james-spec.md` §9).
//!
//! Two encoders, because signatures and public keys want different things.
//!
//! * [`encode_vector`] / [`decode_vector`] carry **signatures**: the whole
//!   vector is one base-`q` integer, so the encoding is exactly
//!   `ceil(n log2 q / 8)` bytes and admits one byte string per value. Without
//!   the canonicity check a signature would be malleable -- over F_2 the
//!   unused high bits of the final byte are free, over odd F_q one may add
//!   `q^n` -- and distinct byte strings would verify against the same message.
//! * [`DigitPacker`] / [`bytes_to_digits`] carry **public keys**, where a
//!   single base conversion over millions of digits would be quadratic.
//!   Digits go in groups of `g` into `B` bytes; the trailing partial group is
//!   packed exactly rather than padded. Decoding validates every group.

use alloc::vec;
use alloc::vec::Vec;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum CodecError {
    /// The byte string is not the length the value requires.
    Length,
    /// The byte string decodes to something outside the value's range, i.e.
    /// it is a second encoding of some other value's bits.
    NonCanonical,
}

// ------------------------------------------------------------ small bignum
//
// Signatures reach 578 base-23 digits, about 2615 bits, so base conversion
// needs more than a machine word. Everything here is O(limbs) per digit and
// runs on public data (a signature is public), so it is kept simple.

pub(crate) fn bn_mul_add(v: &mut Vec<u64>, mul: u64, add: u64) {
    let mut carry = add as u128;
    for w in v.iter_mut() {
        let t = (*w as u128) * (mul as u128) + carry;
        *w = t as u64;
        carry = t >> 64;
    }
    while carry != 0 {
        v.push(carry as u64);
        carry >>= 64;
    }
}

/// `v /= d`, returning the remainder.
pub(crate) fn bn_divmod_small(v: &mut [u64], d: u64) -> u64 {
    let mut rem = 0u128;
    for w in v.iter_mut().rev() {
        let cur = (rem << 64) | (*w as u128);
        *w = (cur / (d as u128)) as u64;
        rem = cur % (d as u128);
    }
    rem as u64
}

pub(crate) fn bn_is_zero(v: &[u64]) -> bool {
    v.iter().all(|&w| w == 0)
}

pub(crate) fn bn_bitlen(v: &[u64]) -> usize {
    for (i, w) in v.iter().enumerate().rev() {
        if *w != 0 {
            return i * 64 + (64 - w.leading_zeros() as usize);
        }
    }
    0
}

/// `q^count` as little-endian limbs.
pub(crate) fn bn_pow(q: u32, count: usize) -> Vec<u64> {
    let mut v = vec![1u64];
    for _ in 0..count {
        bn_mul_add(&mut v, q as u64, 0);
    }
    v
}

/// Bit `i` of a little-endian limb vector.
pub(crate) fn bn_bit(v: &[u64], i: usize) -> bool {
    v.get(i / 64).is_some_and(|w| (w >> (i % 64)) & 1 == 1)
}

/// `a < b` for little-endian limb vectors.
pub(crate) fn bn_lt(a: &[u64], b: &[u64]) -> bool {
    let n = a.len().max(b.len());
    for i in (0..n).rev() {
        let (x, y) = (*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0));
        if x != y {
            return x < y;
        }
    }
    false
}

// ------------------------------------------------------------ vector codec

/// Exact byte length of `count` F_q digits: `ceil(count * log2 q / 8)`.
pub fn vec_bytes(count: usize, q: u32) -> usize {
    if count == 0 {
        return 0;
    }
    if q.is_power_of_two() {
        let k = q.trailing_zeros() as usize;
        return (count * k + 7) / 8;
    }
    // q^count is never a power of two here, so bitlen(q^count - 1) = bitlen(q^count).
    (bn_bitlen(&bn_pow(q, count)) + 7) / 8
}

/// The canonical encoding of a digit vector: one base-`q` integer.
pub fn encode_vector(digits: &[u8], q: u32) -> Vec<u8> {
    let mut v = vec![0u64];
    for &d in digits.iter().rev() {
        bn_mul_add(&mut v, q as u64, d as u64);
    }
    let want = vec_bytes(digits.len(), q);
    let mut out = vec![0u8; want];
    for (i, o) in out.iter_mut().enumerate() {
        *o = (v.get(i / 8).copied().unwrap_or(0) >> (8 * (i % 8))) as u8;
    }
    out
}

/// Inverse of [`encode_vector`], rejecting every non-canonical encoding.
pub fn decode_vector(data: &[u8], count: usize, q: u32) -> Result<Vec<u8>, CodecError> {
    if data.len() != vec_bytes(count, q) {
        return Err(CodecError::Length);
    }
    let mut v = vec![0u64; (data.len() + 7) / 8 + 1];
    for (i, &b) in data.iter().enumerate() {
        v[i / 8] |= (b as u64) << (8 * (i % 8));
    }
    if !bn_lt(&v, &bn_pow(q, count)) {
        return Err(CodecError::NonCanonical);
    }
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        out.push(bn_divmod_small(&mut v, q as u64) as u8);
    }
    debug_assert!(bn_is_zero(&v));
    Ok(out)
}

// ----------------------------------------------------------- grouped codec

/// `(digits per group, bytes per group)`, maximising density.
pub fn pack_shape(q: u32) -> (usize, usize) {
    if q.is_power_of_two() {
        let k = (q.trailing_zeros() as usize).max(1);
        let l = num_lcm(k, 8);
        return (l / k, l / 8);
    }
    let mut best = (0usize, 1usize);
    for b in 1..=16usize {
        let cap: u128 = if b == 16 { u128::MAX } else { 1u128 << (8 * b) };
        let lim = cap / (q as u128);
        let (mut g, mut v) = (0usize, 1u128);
        while v <= lim {
            v *= q as u128;
            g += 1;
        }
        if g > 0 && g * best.1 > best.0 * b {
            best = (g, b);
        }
    }
    best
}

fn num_lcm(a: usize, b: usize) -> usize {
    let (mut x, mut y) = (a, b);
    while y != 0 {
        let t = x % y;
        x = y;
        y = t;
    }
    a / x * b
}

/// Length of the grouped encoding, with the final group packed exactly.
pub fn packed_len(count: usize, q: u32) -> usize {
    let (g, b) = pack_shape(q);
    let (full, rem) = (count / g, count % g);
    full * b + vec_bytes(rem, q)
}

/// Appends F_q digits, yields bytes.
pub struct DigitPacker {
    q: u32,
    g: usize,
    b: usize,
    out: Vec<u8>,
    pending: Vec<u8>,
}

impl DigitPacker {
    pub fn new(q: u32) -> Self {
        let (g, b) = pack_shape(q);
        DigitPacker {
            q,
            g,
            b,
            out: Vec::new(),
            pending: Vec::with_capacity(g),
        }
    }

    pub fn push(&mut self, digits: &[u8]) {
        for &d in digits {
            self.pending.push(d);
            if self.pending.len() == self.g {
                let enc = encode_vector(&self.pending, self.q);
                debug_assert_eq!(enc.len(), self.b);
                self.out.extend_from_slice(&enc);
                self.pending.clear();
            }
        }
    }

    pub fn finish(mut self) -> Vec<u8> {
        if !self.pending.is_empty() {
            let enc = encode_vector(&self.pending, self.q);
            self.out.extend_from_slice(&enc);
        }
        self.out
    }
}

pub fn digits_to_bytes(digits: &[u8], q: u32) -> Vec<u8> {
    let mut p = DigitPacker::new(q);
    p.push(digits);
    p.finish()
}

/// Inverse of [`digits_to_bytes`], rejecting non-canonical encodings.
pub fn bytes_to_digits(data: &[u8], count: usize, q: u32) -> Result<Vec<u8>, CodecError> {
    if data.len() != packed_len(count, q) {
        return Err(CodecError::Length);
    }
    let (g, b) = pack_shape(q);
    let mut out = Vec::with_capacity(count);
    let mut pos = 0usize;
    while out.len() < count {
        let take = core::cmp::min(g, count - out.len());
        let nb = if take == g { b } else { vec_bytes(take, q) };
        let part = decode_vector(&data[pos..pos + nb], take, q)?;
        pos += nb;
        out.extend_from_slice(&part);
    }
    Ok(out)
}

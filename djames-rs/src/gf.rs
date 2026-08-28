//! The base field `F_q` and the extension `K = F_q[t]/(f)`.
//!
//! # Representation
//!
//! An element of `K` is a `Vec<u64>`, packed differently per field so that the
//! public key -- which holds `n(n+1)/2 + n*n_y` of them -- stays in memory:
//!
//! | q | packing | limbs |
//! |---|---|---|
//! | 2 | one coefficient per bit | `ceil(n/64)` |
//! | 4 | two bit-planes, low plane first | `2*ceil(n/64)` |
//! | p odd | one coefficient per byte | `ceil(n/8)` |
//!
//! # Reduction
//!
//! The canonical field polynomials (`d-james-spec.md` §3.2) are sparse with
//! all low exponents below 9, which makes reduction nearly free. Writing
//! `f = t^n + g(t)`, any `H` of degree `>= n` splits as `H = t^n*G + L`, and
//! `t^n = g(t)` gives `H = g*G + L`. Since `deg g < 9`, folding once drops the
//! degree from `2n-2` to at most `n+6`, and a second fold finishes it. The
//! fold count is fixed at construction, so reduction is a straight-line loop
//! with no data-dependent indexing.

use crate::ct;
use alloc::vec;
use alloc::vec::Vec;

pub type Elem = Vec<u64>;

// --------------------------------------------------------------- base field

/// `F_q` for the five `q` the parameter sets use.
///
/// No lookup tables: a 529-entry table indexed by secret values spans several
/// cache lines and would be a timing oracle. Everything is computed.
#[derive(Clone, Copy, Debug)]
pub struct Fq {
    pub q: u32,
    pub p: u32,
    pub k: u32,
    mu: u64,
}

impl Fq {
    pub const fn new(q: u32) -> Self {
        let (p, k) = match q {
            2 => (2, 1),
            4 => (2, 2),
            _ => (q, 1), // 5, 13, 23 are prime
        };
        Fq {
            q,
            p,
            k,
            mu: ct::barrett_mu(q),
        }
    }

    #[inline]
    pub fn add(&self, a: u8, b: u8) -> u8 {
        if self.p == 2 {
            a ^ b // F_2 and F_4 are both characteristic 2, coefficient-wise
        } else {
            let s = (a as u32) + (b as u32);
            let m = ct::lt_mask32(s, self.q);
            ct::select32(m, s, s.wrapping_sub(self.q)) as u8
        }
    }

    #[inline]
    pub fn neg(&self, a: u8) -> u8 {
        if self.p == 2 {
            a
        } else {
            let d = self.q - (a as u32);
            let m = ct::eq_mask32(a as u32, 0);
            ct::select32(m, 0, d) as u8
        }
    }

    #[inline]
    pub fn sub(&self, a: u8, b: u8) -> u8 {
        self.add(a, self.neg(b))
    }

    #[inline]
    pub fn mul(&self, a: u8, b: u8) -> u8 {
        match self.q {
            2 => a & b,
            4 => {
                // F_4 = F_2[u]/(u^2+u+1); a = a0 + a1*u encoded as a0 + 2*a1.
                // (a0+a1u)(b0+b1u) = a0b0 + a1b1 + (a0b1 + a1b0 + a1b1) u
                let (a0, a1) = (a & 1, (a >> 1) & 1);
                let (b0, b1) = (b & 1, (b >> 1) & 1);
                let c0 = (a0 & b0) ^ (a1 & b1);
                let c1 = (a0 & b1) ^ (a1 & b0) ^ (a1 & b1);
                c0 | (c1 << 1)
            }
            _ => ct::barrett_u32((a as u32) * (b as u32), self.q, self.mu) as u8,
        }
    }

    /// Multiplicative inverse; `inv(0)` is 0, which callers must not rely on.
    pub fn inv(&self, a: u8) -> u8 {
        // a^(q-2) by square-and-multiply over a fixed exponent: the schedule
        // depends only on q, never on a.
        let mut r = 1u8;
        let mut b = a;
        let mut e = self.q - 2;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            b = self.mul(b, b);
            e >>= 1;
        }
        r
    }
}

// ------------------------------------------------------------- clmul (q = 2)

/// Carry-less product of two 32-bit values, constant time.
///
/// The four-way nibble split keeps every output nibble below 16, so an
/// ordinary integer multiply never carries between the classes we keep. This
/// is BearSSL's portable GHASH technique.
#[inline]
fn bmul32(x: u32, y: u32) -> u64 {
    const M0: u64 = 0x1111_1111_1111_1111;
    const M1: u64 = 0x2222_2222_2222_2222;
    const M2: u64 = 0x4444_4444_4444_4444;
    const M3: u64 = 0x8888_8888_8888_8888;
    let (x, y) = (x as u64, y as u64);
    let (x0, x1, x2, x3) = (x & M0, x & M1, x & M2, x & M3);
    let (y0, y1, y2, y3) = (y & M0, y & M1, y & M2, y & M3);
    let z0 = x0.wrapping_mul(y0) ^ x1.wrapping_mul(y3) ^ x2.wrapping_mul(y2) ^ x3.wrapping_mul(y1);
    let z1 = x0.wrapping_mul(y1) ^ x1.wrapping_mul(y0) ^ x2.wrapping_mul(y3) ^ x3.wrapping_mul(y2);
    let z2 = x0.wrapping_mul(y2) ^ x1.wrapping_mul(y1) ^ x2.wrapping_mul(y0) ^ x3.wrapping_mul(y3);
    let z3 = x0.wrapping_mul(y3) ^ x1.wrapping_mul(y2) ^ x2.wrapping_mul(y1) ^ x3.wrapping_mul(y0);
    (z0 & M0) | (z1 & M1) | (z2 & M2) | (z3 & M3)
}

/// Carry-less product of two 64-bit values -> 128 bits, via Karatsuba.
#[inline]
pub fn clmul64(x: u64, y: u64) -> (u64, u64) {
    let (a0, a1) = (x as u32, (x >> 32) as u32);
    let (b0, b1) = (y as u32, (y >> 32) as u32);
    let z0 = bmul32(a0, b0);
    let z2 = bmul32(a1, b1);
    let z1 = bmul32(a0 ^ a1, b0 ^ b1) ^ z0 ^ z2;
    (z0 ^ (z1 << 32), z2 ^ (z1 >> 32))
}

// ---------------------------------------------------------------- the field

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    /// One coefficient per bit.
    Bits,
    /// Two bit-planes.
    Planes,
    /// One coefficient per byte.
    Bytes,
}

/// `K = F_q[t] / (f)`, with `alpha` the class of `t`.
pub struct Ext {
    pub q: u32,
    pub n: usize,
    pub fq: Fq,
    kind: Kind,
    /// u64 words in one element.
    pub limbs: usize,
    /// words in one *plane* (Planes kind) or in the whole element otherwise.
    plane: usize,
    /// `f = t^n + sum c_i t^i`, low coefficients.
    #[allow(dead_code)]
    fpoly: Vec<u8>,
    /// nonzero `(i, c_i)` of the low part.
    fsupport: Vec<(usize, u8)>,
    /// how many folds reduce degree `2n-2` below `n`.
    folds: usize,
    /// Frobenius tables: `(k, [(alpha^i)^(q^k)]_{i<n})`.
    frob: Vec<(usize, Vec<Elem>)>,
}

impl Ext {
    pub fn new(q: u32, n: usize, fpoly: &[u8]) -> Self {
        assert_eq!(fpoly.len(), n);
        let fq = Fq::new(q);
        let (kind, plane, limbs) = match q {
            2 => (Kind::Bits, (n + 63) / 64, (n + 63) / 64),
            4 => (Kind::Planes, (n + 63) / 64, 2 * ((n + 63) / 64)),
            _ => (Kind::Bytes, (n + 7) / 8, (n + 7) / 8),
        };
        let fsupport: Vec<(usize, u8)> = fpoly
            .iter()
            .enumerate()
            .filter(|(_, &c)| c != 0)
            .map(|(i, &c)| (i, c))
            .collect();
        let top = fsupport.last().map(|&(i, _)| i).unwrap_or(0);
        // One fold takes degree 2n-2 to at most (n-2)+top; repeat until < n.
        let mut folds = 0usize;
        let mut deg = 2 * n - 2;
        while deg >= n {
            deg = (deg - n) + top;
            folds += 1;
            assert!(folds < 8, "field polynomial is not sparse enough");
        }
        let mut e = Ext {
            q,
            n,
            fq,
            kind,
            limbs,
            plane,
            fpoly: fpoly.to_vec(),
            fsupport,
            folds,
            frob: Vec::new(),
        };
        e.build_frobenius();
        e
    }

    // -- element construction -------------------------------------------
    pub fn zero(&self) -> Elem {
        vec![0u64; self.limbs]
    }

    pub fn one(&self) -> Elem {
        let mut e = self.zero();
        self.set_coef(&mut e, 0, 1);
        e
    }

    pub fn alpha(&self) -> Elem {
        let mut e = self.zero();
        self.set_coef(&mut e, 1, 1);
        e
    }

    #[inline]
    pub fn coef(&self, a: &[u64], i: usize) -> u8 {
        match self.kind {
            Kind::Bits => ((a[i / 64] >> (i % 64)) & 1) as u8,
            Kind::Planes => {
                let lo = ((a[i / 64] >> (i % 64)) & 1) as u8;
                let hi = ((a[self.plane + i / 64] >> (i % 64)) & 1) as u8;
                lo | (hi << 1)
            }
            Kind::Bytes => ((a[i / 8] >> (8 * (i % 8))) & 0xff) as u8,
        }
    }

    #[inline]
    pub fn set_coef(&self, a: &mut [u64], i: usize, v: u8) {
        match self.kind {
            Kind::Bits => {
                let m = 1u64 << (i % 64);
                a[i / 64] = (a[i / 64] & !m) | (((v & 1) as u64) << (i % 64));
            }
            Kind::Planes => {
                let m = 1u64 << (i % 64);
                a[i / 64] = (a[i / 64] & !m) | (((v & 1) as u64) << (i % 64));
                a[self.plane + i / 64] =
                    (a[self.plane + i / 64] & !m) | ((((v >> 1) & 1) as u64) << (i % 64));
            }
            Kind::Bytes => {
                let sh = 8 * (i % 8);
                a[i / 8] = (a[i / 8] & !(0xffu64 << sh)) | ((v as u64) << sh);
            }
        }
    }

    pub fn coords(&self, a: &[u64]) -> Vec<u8> {
        (0..self.n).map(|i| self.coef(a, i)).collect()
    }

    pub fn from_coords(&self, cs: &[u8]) -> Elem {
        let mut e = self.zero();
        for (i, &c) in cs.iter().enumerate().take(self.n) {
            self.set_coef(&mut e, i, c);
        }
        e
    }

    pub fn is_zero(&self, a: &[u64]) -> bool {
        a.iter().all(|&w| w == 0)
    }

    // -- additive structure ----------------------------------------------
    pub fn add(&self, a: &[u64], b: &[u64]) -> Elem {
        let mut r = a.to_vec();
        self.add_assign(&mut r, b);
        r
    }

    pub fn add_assign(&self, a: &mut [u64], b: &[u64]) {
        match self.kind {
            Kind::Bits | Kind::Planes => {
                for (x, &y) in a.iter_mut().zip(b.iter()) {
                    *x ^= y;
                }
            }
            Kind::Bytes => {
                for i in 0..self.n {
                    let v = self.fq.add(self.coef(a, i), self.coef(b, i));
                    self.set_coef(a, i, v);
                }
            }
        }
    }

    pub fn neg(&self, a: &[u64]) -> Elem {
        if self.fq.p == 2 {
            return a.to_vec();
        }
        let mut r = self.zero();
        for i in 0..self.n {
            self.set_coef(&mut r, i, self.fq.neg(self.coef(a, i)));
        }
        r
    }

    pub fn sub(&self, a: &[u64], b: &[u64]) -> Elem {
        if self.fq.p == 2 {
            return self.add(a, b);
        }
        self.add(a, &self.neg(b))
    }

    /// Multiply by a scalar of `F_q`.
    pub fn scal(&self, a: &[u64], s: u8) -> Elem {
        match self.kind {
            Kind::Bits => {
                let m = ct::mask64((s & 1) as u64);
                a.iter().map(|&w| w & m).collect()
            }
            _ => {
                let mut r = self.zero();
                for i in 0..self.n {
                    self.set_coef(&mut r, i, self.fq.mul(self.coef(a, i), s));
                }
                r
            }
        }
    }

    pub fn eq(&self, a: &[u64], b: &[u64]) -> bool {
        a == b
    }
}

// ------------------------------------------------------- bit-level helpers

/// `dst ^= src << bits`, over packed bit words.
fn xor_shl(dst: &mut [u64], src: &[u64], bits: usize) {
    let (w, b) = (bits / 64, bits % 64);
    if b == 0 {
        for i in 0..src.len() {
            if i + w < dst.len() {
                dst[i + w] ^= src[i];
            }
        }
    } else {
        for i in 0..src.len() {
            if i + w < dst.len() {
                dst[i + w] ^= src[i] << b;
            }
            if i + w + 1 < dst.len() {
                dst[i + w + 1] ^= src[i] >> (64 - b);
            }
        }
    }
}

/// `out = src >> bits`, over packed bit words.
fn shr_into(out: &mut [u64], src: &[u64], bits: usize) {
    let (w, b) = (bits / 64, bits % 64);
    for o in out.iter_mut() {
        *o = 0;
    }
    for i in 0..out.len() {
        if i + w < src.len() {
            out[i] = src[i + w] >> b;
            if b != 0 && i + w + 1 < src.len() {
                out[i] |= src[i + w + 1] << (64 - b);
            }
        }
    }
}

/// Clear every bit at index `>= n`.
fn mask_below(p: &mut [u64], n: usize) {
    let (w, b) = (n / 64, n % 64);
    if b != 0 && w < p.len() {
        p[w] &= (1u64 << b) - 1;
    }
    let first = if b == 0 { w } else { w + 1 };
    for x in p.iter_mut().skip(first) {
        *x = 0;
    }
}

impl Ext {
    // -- multiplication --------------------------------------------------

    pub fn mul(&self, a: &[u64], b: &[u64]) -> Elem {
        match self.kind {
            Kind::Bits => self.mul_bits(a, b),
            Kind::Planes => self.mul_planes(a, b),
            Kind::Bytes => self.mul_bytes(a, b),
        }
    }

    /// Carry-less polynomial product of two bit-packed operands.
    fn clmul_poly(&self, a: &[u64], b: &[u64], out: &mut [u64]) {
        let l = self.plane;
        for o in out.iter_mut() {
            *o = 0;
        }
        for i in 0..l {
            for j in 0..l {
                let (lo, hi) = clmul64(a[i], b[j]);
                out[i + j] ^= lo;
                out[i + j + 1] ^= hi;
            }
        }
    }

    /// Fold a bit-packed polynomial of degree `< 2n-1` back below `t^n`.
    ///
    /// `f = t^n + g`, so `H = t^n*G + L` reduces to `g*G + L`. `folds` is
    /// fixed at construction, so this loop is straight-line.
    fn reduce_bits_in(&self, p: &mut [u64]) {
        let mut g = vec![0u64; p.len()];
        for _ in 0..self.folds {
            shr_into(&mut g, p, self.n);
            mask_below(p, self.n);
            for &(i, _) in &self.fsupport {
                xor_shl(p, &g, i); // over F_2 every c_i is 1
            }
        }
    }

    fn mul_bits(&self, a: &[u64], b: &[u64]) -> Elem {
        let mut p = vec![0u64; 2 * self.plane + 1];
        self.clmul_poly(a, b, &mut p);
        self.reduce_bits_in(&mut p);
        p.truncate(self.limbs);
        p
    }

    /// `F_4[t]` product: three `F_2[t]` products via Karatsuba, recombined
    /// through `u^2 = u + 1`.
    fn mul_planes(&self, a: &[u64], b: &[u64]) -> Elem {
        let l = self.plane;
        let wide = 2 * l + 1;
        let (a0, a1) = (&a[..l], &a[l..2 * l]);
        let (b0, b1) = (&b[..l], &b[l..2 * l]);
        let ax: Vec<u64> = a0.iter().zip(a1).map(|(x, y)| x ^ y).collect();
        let bx: Vec<u64> = b0.iter().zip(b1).map(|(x, y)| x ^ y).collect();
        let (mut t0, mut t1, mut t2) = (vec![0u64; wide], vec![0u64; wide], vec![0u64; wide]);
        self.clmul_poly(a0, b0, &mut t0);
        self.clmul_poly(a1, b1, &mut t1);
        self.clmul_poly(&ax, &bx, &mut t2);
        // c0 = t0 + t1 ; c1 = (a0b1 + a1b0) + t1 = (t2 + t0 + t1) + t1 = t2 + t0
        let mut c0: Vec<u64> = t0.iter().zip(&t1).map(|(x, y)| x ^ y).collect();
        let mut c1: Vec<u64> = t2.iter().zip(&t0).map(|(x, y)| x ^ y).collect();
        self.reduce_planes_in(&mut c0, &mut c1);
        let mut r = vec![0u64; self.limbs];
        r[..l].copy_from_slice(&c0[..l]);
        r[l..2 * l].copy_from_slice(&c1[..l]);
        r
    }

    fn reduce_planes_in(&self, p0: &mut [u64], p1: &mut [u64]) {
        let mut g0 = vec![0u64; p0.len()];
        let mut g1 = vec![0u64; p1.len()];
        for _ in 0..self.folds {
            shr_into(&mut g0, p0, self.n);
            shr_into(&mut g1, p1, self.n);
            mask_below(p0, self.n);
            mask_below(p1, self.n);
            for &(i, c) in &self.fsupport {
                // scale (g0, g1) by c in F_4, then shift into place
                let (s0, s1): (Vec<u64>, Vec<u64>) = match c {
                    1 => (g0.clone(), g1.clone()),
                    2 => (g1.clone(), g0.iter().zip(&g1).map(|(x, y)| x ^ y).collect()),
                    _ => (g0.iter().zip(&g1).map(|(x, y)| x ^ y).collect(), g0.clone()),
                };
                xor_shl(p0, &s0, i);
                xor_shl(p1, &s1, i);
            }
        }
    }

    /// Schoolbook product over odd `F_p`.
    ///
    /// Coefficients are unpacked to `u32`, accumulated without reduction --
    /// the largest partial sum is `n*(q-1)^2`, far below 2^32 -- then reduced
    /// once and folded modulo `f` from the top down. No coefficient test is
    /// skipped, so the work does not depend on the values.
    fn mul_bytes(&self, a: &[u64], b: &[u64]) -> Elem {
        let n = self.n;
        let (ac, bc) = (self.coords(a), self.coords(b));
        let mut acc = vec![0u32; 2 * n];
        for i in 0..n {
            let ai = ac[i] as u32;
            let dst = &mut acc[i..i + n];
            for (d, &bj) in dst.iter_mut().zip(bc.iter()) {
                *d += ai * (bj as u32);
            }
        }
        let mut c = vec![0u8; 2 * n];
        let mu = ct::barrett_mu(self.q);
        for (o, &v) in c.iter_mut().zip(acc.iter()) {
            *o = ct::barrett_u32(v, self.q, mu) as u8;
        }
        // Fold from the top: every target index is below the source, so one
        // downward pass suffices.
        for k in (n..2 * n - 1).rev() {
            let v = c[k];
            c[k] = 0;
            for &(i, ci) in &self.fsupport {
                let t = self.fq.mul(v, ci);
                let dst = k - n + i;
                c[dst] = self.fq.sub(c[dst], t); // t^n = -(sum c_i t^i)
            }
        }
        self.from_coords(&c[..n])
    }

    // -- exponentiation and inversion ------------------------------------

    /// `a^e` for a *public* exponent.
    pub fn pow(&self, a: &[u64], e: u64) -> Elem {
        let mut r = self.one();
        let mut b = a.to_vec();
        let mut e = e;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(&r, &b);
            }
            e >>= 1;
            if e > 0 {
                b = self.mul(&b, &b);
            }
        }
        r
    }

    /// Multiplicative inverse.
    ///
    /// Via the norm rather than `a^(q^n - 2)`, whose exponent does not fit in
    /// a machine word: with `B = a^(q + q^2 + ... + q^(n-1))`, the product
    /// `a*B` is the norm and lies in `F_q`, so `a^-1 = B * norm^-1`.
    pub fn inv(&self, a: &[u64]) -> Elem {
        debug_assert!(!self.is_zero(a));
        let aq = self.pow(a, self.q as u64);
        let mut bb = aq.clone();
        for _ in 1..self.n - 1 {
            bb = self.pow(&bb, self.q as u64);
            bb = self.mul(&bb, &aq);
        }
        let norm = self.mul(a, &bb);
        let s = self.fq.inv(self.coef(&norm, 0));
        self.scal(&bb, s)
    }

    // -- Frobenius --------------------------------------------------------

    fn build_frobenius(&mut self) {
        // Exponents needed: 0..=8 covers every `d` in the parameter sets, and
        // the binary-powering ladder in poly::x_pow_qn visits exactly the
        // prefixes of n written in binary.
        let mut ks: Vec<usize> = (0..=8).collect();
        let mut e = 0usize;
        for bit in (0..usize::BITS - self.n.leading_zeros()).rev() {
            e = e * 2 + ((self.n >> bit) & 1);
            ks.push(e);
            ks.push(e * 2);
        }
        ks.sort_unstable();
        ks.dedup();
        ks.retain(|&k| k <= self.n);
        for k in ks {
            let mut beta = self.alpha();
            for _ in 0..k {
                beta = self.pow(&beta, self.q as u64);
            }
            let mut tbl = Vec::with_capacity(self.n);
            let mut cur = self.one();
            for _ in 0..self.n {
                tbl.push(cur.clone());
                cur = self.mul(&cur, &beta);
            }
            self.frob.push((k, tbl));
        }
    }

    /// `a^(q^k)`.
    pub fn frob(&self, a: &[u64], k: usize) -> Elem {
        let k = k % self.n;
        if k == 0 {
            return a.to_vec();
        }
        let tbl = &self
            .frob
            .iter()
            .find(|(kk, _)| *kk == k)
            .unwrap_or_else(|| panic!("no Frobenius table for q^{k}"))
            .1;
        let mut r = self.zero();
        for i in 0..self.n {
            let c = self.coef(a, i);
            if self.kind == Kind::Bits {
                let m = ct::mask64(c as u64);
                for (x, &y) in r.iter_mut().zip(tbl[i].iter()) {
                    *x ^= y & m;
                }
            } else {
                let t = self.scal(&tbl[i], c);
                self.add_assign(&mut r, &t);
            }
        }
        r
    }
}

// ------------------------------------------------- packed F_q^len vectors
//
// The public key holds `n(n+1)/2 + n*n_y` coefficient vectors of length m --
// 167331 of them at `james-256-q2`. One byte per coefficient would need 86 MB;
// the packing below needs 11 MB.

/// A packed vector over `F_q` of fixed length, with no multiplication.
pub struct PackedVec {
    pub q: u32,
    pub len: usize,
    pub fq: Fq,
    kind: Kind,
    plane: usize,
    pub words: usize,
}

impl PackedVec {
    pub fn new(q: u32, len: usize) -> Self {
        let fq = Fq::new(q);
        let (kind, plane, words) = match q {
            2 => (Kind::Bits, (len + 63) / 64, (len + 63) / 64),
            4 => (Kind::Planes, (len + 63) / 64, 2 * ((len + 63) / 64)),
            _ => (Kind::Bytes, (len + 7) / 8, (len + 7) / 8),
        };
        PackedVec {
            q,
            len,
            fq,
            kind,
            plane,
            words,
        }
    }

    pub fn zero(&self) -> Vec<u64> {
        vec![0u64; self.words]
    }

    #[inline]
    pub fn coef(&self, a: &[u64], i: usize) -> u8 {
        match self.kind {
            Kind::Bits => ((a[i / 64] >> (i % 64)) & 1) as u8,
            Kind::Planes => {
                let lo = ((a[i / 64] >> (i % 64)) & 1) as u8;
                let hi = ((a[self.plane + i / 64] >> (i % 64)) & 1) as u8;
                lo | (hi << 1)
            }
            Kind::Bytes => ((a[i / 8] >> (8 * (i % 8))) & 0xff) as u8,
        }
    }

    #[inline]
    pub fn set_coef(&self, a: &mut [u64], i: usize, v: u8) {
        match self.kind {
            Kind::Bits => {
                let m = 1u64 << (i % 64);
                a[i / 64] = (a[i / 64] & !m) | (((v & 1) as u64) << (i % 64));
            }
            Kind::Planes => {
                let m = 1u64 << (i % 64);
                a[i / 64] = (a[i / 64] & !m) | (((v & 1) as u64) << (i % 64));
                a[self.plane + i / 64] =
                    (a[self.plane + i / 64] & !m) | ((((v >> 1) & 1) as u64) << (i % 64));
            }
            Kind::Bytes => {
                let sh = 8 * (i % 8);
                a[i / 8] = (a[i / 8] & !(0xffu64 << sh)) | ((v as u64) << sh);
            }
        }
    }

    pub fn from_digits(&self, cs: &[u8]) -> Vec<u64> {
        let mut v = self.zero();
        for (i, &c) in cs.iter().enumerate().take(self.len) {
            self.set_coef(&mut v, i, c);
        }
        v
    }

    pub fn to_digits(&self, a: &[u64]) -> Vec<u8> {
        (0..self.len).map(|i| self.coef(a, i)).collect()
    }

    pub fn is_zero(&self, a: &[u64]) -> bool {
        a.iter().all(|&w| w == 0)
    }

    /// `acc += s * v`.
    pub fn scal_add_assign(&self, acc: &mut [u64], v: &[u64], s: u8) {
        match self.kind {
            Kind::Bits => {
                let m = ct::mask64((s & 1) as u64);
                for (x, &y) in acc.iter_mut().zip(v.iter()) {
                    *x ^= y & m;
                }
            }
            _ => {
                for i in 0..self.len {
                    let t = self.fq.mul(self.coef(v, i), s);
                    let cur = self.coef(acc, i);
                    self.set_coef(acc, i, self.fq.add(cur, t));
                }
            }
        }
    }
}

// ---------------------------------------------------- allocation-free path
//
// Public-key assembly performs O(n^2) products per term. Allocating a fresh
// element for each one dominates; these variants write through caller-owned
// buffers instead.

/// Reusable working buffers for [`Ext::mul_into`].
pub struct Scratch {
    wide: Vec<u64>,
    w2: Vec<u64>,
    w3: Vec<u64>,
    g0: Vec<u64>,
    g1: Vec<u64>,
    acc: Vec<u32>,
    ca: Vec<u8>,
    cb: Vec<u8>,
    cc: Vec<u8>,
}

impl Ext {
    pub fn scratch(&self) -> Scratch {
        let wide = 2 * self.plane + 1;
        Scratch {
            wide: vec![0u64; wide],
            w2: vec![0u64; wide],
            w3: vec![0u64; wide],
            g0: vec![0u64; wide],
            g1: vec![0u64; wide],
            acc: vec![0u32; 2 * self.n],
            ca: vec![0u8; self.n],
            cb: vec![0u8; self.n],
            cc: vec![0u8; 2 * self.n],
        }
    }

    /// `out = a * b`, allocating nothing.
    pub fn mul_into(&self, out: &mut [u64], a: &[u64], b: &[u64], s: &mut Scratch) {
        match self.kind {
            Kind::Bits => {
                self.clmul_poly(a, b, &mut s.wide);
                self.reduce_bits_scratch(&mut s.wide, &mut s.g0);
                out.copy_from_slice(&s.wide[..self.limbs]);
            }
            Kind::Planes => {
                let l = self.plane;
                let (a0, a1) = (&a[..l], &a[l..2 * l]);
                let (b0, b1) = (&b[..l], &b[l..2 * l]);
                for i in 0..l {
                    s.g0[i] = a0[i] ^ a1[i];
                    s.g1[i] = b0[i] ^ b1[i];
                }
                let (ax, bx) = (s.g0[..l].to_vec(), s.g1[..l].to_vec());
                self.clmul_poly(a0, b0, &mut s.wide); // t0
                self.clmul_poly(a1, b1, &mut s.w2); // t1
                self.clmul_poly(&ax, &bx, &mut s.w3); // t2
                for i in 0..s.wide.len() {
                    let (t0, t1, t2) = (s.wide[i], s.w2[i], s.w3[i]);
                    s.wide[i] = t0 ^ t1; // c0
                    s.w3[i] = t2 ^ t0; // c1
                }
                let (mut c0, mut c1) = (core::mem::take(&mut s.wide), core::mem::take(&mut s.w3));
                self.reduce_planes_in(&mut c0, &mut c1);
                out[..l].copy_from_slice(&c0[..l]);
                out[l..2 * l].copy_from_slice(&c1[..l]);
                s.wide = c0;
                s.w3 = c1;
            }
            Kind::Bytes => {
                let n = self.n;
                for i in 0..n {
                    s.ca[i] = self.coef(a, i);
                    s.cb[i] = self.coef(b, i);
                }
                for v in s.acc.iter_mut() {
                    *v = 0;
                }
                for i in 0..n {
                    let ai = s.ca[i] as u32;
                    let dst = &mut s.acc[i..i + n];
                    for (d, &bj) in dst.iter_mut().zip(s.cb.iter()) {
                        *d += ai * (bj as u32);
                    }
                }
                let mu = ct::barrett_mu(self.q);
                for (o, &v) in s.cc.iter_mut().zip(s.acc.iter()) {
                    *o = ct::barrett_u32(v, self.q, mu) as u8;
                }
                for k in (n..2 * n - 1).rev() {
                    let v = s.cc[k];
                    s.cc[k] = 0;
                    for &(i, ci) in &self.fsupport {
                        let t = self.fq.mul(v, ci);
                        let dst = k - n + i;
                        s.cc[dst] = self.fq.sub(s.cc[dst], t);
                    }
                }
                for w in out.iter_mut() {
                    *w = 0;
                }
                for i in 0..n {
                    self.set_coef(out, i, s.cc[i]);
                }
            }
        }
    }

    fn reduce_bits_scratch(&self, p: &mut [u64], g: &mut [u64]) {
        for _ in 0..self.folds {
            shr_into(g, p, self.n);
            mask_below(p, self.n);
            for &(i, _) in &self.fsupport {
                xor_shl(p, g, i);
            }
        }
    }

    /// Project through an `F_q`-linear map given by its `n` packed columns:
    /// `out += sum_c coef(e, c) * cols[c]`.
    ///
    /// Specialised to avoid materialising the coordinate vector; the F_2 path
    /// is a masked XOR, so it stays branch-free.
    pub fn project_into(&self, out: &mut [u64], e: &[u64], cols: &[Vec<u64>], w: &PackedVec) {
        match self.kind {
            Kind::Bits => {
                for c in 0..self.n {
                    let m = ct::mask64((e[c / 64] >> (c % 64)) & 1);
                    for (x, &y) in out.iter_mut().zip(cols[c].iter()) {
                        *x ^= y & m;
                    }
                }
            }
            _ => {
                for c in 0..self.n {
                    w.scal_add_assign(out, &cols[c], self.coef(e, c));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params;

    fn ext(q: u32, n: usize) -> Ext {
        let p = params::ALL.iter().find(|p| p.q == q && p.n == n).unwrap();
        Ext::new(q, n, p.fpoly)
    }

    #[test]
    fn base_field_axioms() {
        for q in [2u32, 4, 5, 13, 23] {
            let f = Fq::new(q);
            for a in 0..q as u8 {
                assert_eq!(f.add(a, f.neg(a)), 0, "q={q} a={a}");
                if a != 0 {
                    assert_eq!(f.mul(a, f.inv(a)), 1, "q={q} a={a}");
                }
                for b in 0..q as u8 {
                    for c in 0..q as u8 {
                        assert_eq!(
                            f.mul(a, f.add(b, c)),
                            f.add(f.mul(a, b), f.mul(a, c)),
                            "q={q}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn clmul_matches_reference() {
        // Bit-serial carry-less multiply, obviously correct, as the oracle.
        fn slow(x: u64, y: u64) -> (u64, u64) {
            let (mut lo, mut hi) = (0u64, 0u64);
            for i in 0..64 {
                if (y >> i) & 1 == 1 {
                    lo ^= x << i;
                    if i > 0 {
                        hi ^= x >> (64 - i);
                    }
                }
            }
            (lo, hi)
        }
        let mut s = 0x1234_5678_9abc_def0u64;
        for _ in 0..2000 {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let t = s.rotate_left(17) ^ 0xa5a5_5a5a_1234_4321;
            assert_eq!(clmul64(s, t), slow(s, t));
        }
        assert_eq!(clmul64(0, 12345), (0, 0));
        assert_eq!(clmul64(1, 12345), (12345, 0));
    }

    /// Cross-check against `djames-py`: same field polynomial, same inputs.
    #[test]
    fn matches_python_reference() {
        /// `(q, n, a, b, a*b, a^-1, a^(q^3))`, all from `djames-py`.
        type Case = (
            u32,
            usize,
            &'static [u8],
            &'static [u8],
            &'static [u8],
            &'static [u8],
            &'static [u8],
        );
        let cases: &[Case] = &[
            (
                2,
                32,
                &[
                    1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
                    1, 0, 1, 0, 1, 0,
                ],
                &[
                    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                    0, 1, 0, 1, 0, 1,
                ],
                &[
                    0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
                    0, 1, 0, 0, 0, 1,
                ],
                &[
                    1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
                    0, 1, 1, 1, 0, 0,
                ],
                &[
                    1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1,
                    0, 0, 0, 1, 1, 1,
                ],
            ),
            (
                4,
                24,
                &[
                    1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2,
                ],
                &[
                    2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1,
                ],
                &[
                    1, 0, 3, 0, 0, 1, 1, 0, 1, 0, 2, 2, 0, 2, 3, 0, 1, 0, 2, 2, 0, 2, 3, 0,
                ],
                &[
                    3, 0, 1, 2, 0, 0, 0, 3, 1, 3, 2, 2, 2, 3, 1, 2, 3, 1, 1, 0, 1, 1, 3, 2,
                ],
                &[
                    0, 1, 3, 3, 3, 0, 0, 1, 0, 0, 2, 2, 3, 1, 1, 1, 0, 2, 1, 3, 2, 3, 2, 3,
                ],
            ),
            (
                5,
                21,
                &[
                    1, 4, 2, 0, 3, 1, 4, 2, 0, 3, 1, 4, 2, 0, 3, 1, 4, 2, 0, 3, 1,
                ],
                &[
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                ],
                &[
                    2, 3, 3, 4, 1, 4, 3, 3, 4, 1, 4, 3, 3, 4, 1, 4, 3, 3, 4, 1, 4,
                ],
                &[
                    4, 4, 3, 1, 2, 4, 3, 1, 2, 4, 3, 1, 2, 4, 3, 1, 2, 4, 3, 1, 4,
                ],
                &[
                    1, 4, 1, 3, 1, 3, 0, 3, 3, 0, 1, 3, 3, 3, 3, 0, 3, 3, 3, 3, 0,
                ],
            ),
            (
                13,
                18,
                &[1, 4, 7, 10, 0, 3, 6, 9, 12, 2, 5, 8, 11, 1, 4, 7, 10, 0],
                &[2, 7, 12, 4, 9, 1, 6, 11, 3, 8, 0, 5, 10, 2, 7, 12, 4, 9],
                &[11, 4, 4, 4, 10, 2, 12, 7, 6, 2, 1, 9, 6, 11, 4, 4, 4, 10],
                &[1, 5, 8, 7, 11, 8, 7, 11, 8, 7, 11, 8, 7, 11, 8, 7, 11, 8],
                &[1, 3, 8, 3, 0, 4, 6, 10, 10, 11, 6, 2, 11, 4, 12, 6, 12, 0],
            ),
            (
                23,
                16,
                &[1, 4, 7, 10, 13, 16, 19, 22, 2, 5, 8, 11, 14, 17, 20, 0],
                &[2, 7, 12, 17, 22, 4, 9, 14, 19, 1, 6, 11, 16, 21, 3, 8],
                &[14, 18, 7, 12, 18, 10, 19, 7, 5, 21, 17, 1, 4, 11, 7, 0],
                &[6, 1, 12, 19, 9, 20, 1, 15, 18, 17, 2, 7, 13, 11, 4, 14],
                &[9, 8, 22, 11, 14, 20, 17, 5, 22, 13, 14, 13, 0, 17, 22, 7],
            ),
        ];
        for &(q, n, av, bv, prod, inv, fr) in cases {
            let k = ext(q, n);
            let a = k.from_coords(av);
            let b = k.from_coords(bv);
            assert_eq!(k.coords(&k.mul(&a, &b)), prod, "mul q={q}");
            assert_eq!(k.coords(&k.inv(&a)), inv, "inv q={q}");
            assert_eq!(k.coords(&k.frob(&a, 3)), fr, "frob q={q}");
        }
    }

    #[test]
    fn field_axioms_over_k() {
        for (q, n) in [(2u32, 32usize), (4, 24), (5, 21), (13, 18), (23, 16)] {
            let k = ext(q, n);
            let mut s = 12345u64;
            let mut rnd = || {
                let cs: Vec<u8> = (0..n)
                    .map(|_| {
                        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((s >> 33) % q as u64) as u8
                    })
                    .collect();
                k.from_coords(&cs)
            };
            for _ in 0..20 {
                let (a, b, c) = (rnd(), rnd(), rnd());
                assert!(k.eq(&k.mul(&a, &b), &k.mul(&b, &a)));
                assert!(k.eq(&k.mul(&a, &k.mul(&b, &c)), &k.mul(&k.mul(&a, &b), &c)));
                assert!(k.eq(
                    &k.mul(&a, &k.add(&b, &c)),
                    &k.add(&k.mul(&a, &b), &k.mul(&a, &c))
                ));
                assert!(k.eq(&k.sub(&k.add(&a, &b), &b), &a));
                if !k.is_zero(&a) {
                    assert!(k.eq(&k.mul(&a, &k.inv(&a)), &k.one()));
                }
                // a^(q^n) = a, and frob composes
                assert!(k.eq(&k.frob(&a, n), &a));
                let mut r = a.clone();
                for _ in 0..3 {
                    r = k.pow(&r, q as u64);
                }
                assert!(k.eq(&k.frob(&a, 3), &r));
            }
        }
    }
}

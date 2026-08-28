//! Univariate polynomials over `K = F_q^n`, and root finding.
//!
//! Inverting the central map means finding a root in `K` of a polynomial of
//! degree `D = q^d + 1`, which the parameter sets keep at 24 or below. So
//! everything here optimises for "tiny degree, huge coefficient field":
//! coefficient vectors, schoolbook products, and a Frobenius ladder that
//! never exponentiates by more than `q`.
//!
//! A polynomial is a `Vec<Elem>` indexed by degree with no trailing zero; the
//! zero polynomial is empty.
//!
//! # Constant time
//!
//! **This module is not constant time, and cannot cheaply be made so.** The
//! number of roots and the control flow of equal-degree splitting depend on
//! secret-derived data. That is inherent to HFE inversion rather than
//! particular to this code. The arithmetic underneath it (`gf`) is constant
//! time, so what leaks is the shape of the factorisation, not the field
//! operations.

use crate::codec::{bn_bit, bn_bitlen, bn_divmod_small, bn_pow, bn_sub1};
use crate::gf::{Elem, Ext};
use crate::symmetric::{sample_fq, Xof};
use alloc::vec;
use alloc::vec::Vec;

pub type Poly = Vec<Elem>;

pub fn norm(k: &Ext, mut p: Poly) -> Poly {
    while let Some(last) = p.last() {
        if k.is_zero(last) {
            p.pop();
        } else {
            break;
        }
    }
    p
}

/// Degree, or `-1` for the zero polynomial.
pub fn deg(p: &Poly) -> isize {
    p.len() as isize - 1
}

pub fn add(k: &Ext, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let x = a.get(i).cloned().unwrap_or_else(|| k.zero());
        match b.get(i) {
            Some(y) => out.push(k.add(&x, y)),
            None => out.push(x),
        }
    }
    norm(k, out)
}

pub fn sub(k: &Ext, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let x = a.get(i).cloned().unwrap_or_else(|| k.zero());
        match b.get(i) {
            Some(y) => out.push(k.sub(&x, y)),
            None => out.push(x),
        }
    }
    norm(k, out)
}

pub fn scal(k: &Ext, a: &Poly, s: &Elem) -> Poly {
    if k.is_zero(s) {
        return Vec::new();
    }
    norm(k, a.iter().map(|c| k.mul(c, s)).collect())
}

pub fn mul(k: &Ext, a: &Poly, b: &Poly) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![k.zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if k.is_zero(x) {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if !k.is_zero(y) {
                let t = k.mul(x, y);
                out[i + j] = k.add(&out[i + j], &t);
            }
        }
    }
    norm(k, out)
}

pub fn monic(k: &Ext, a: &Poly) -> Poly {
    match a.last() {
        None => Vec::new(),
        Some(lc) if k.eq(lc, &k.one()) => a.clone(),
        Some(lc) => {
            let inv = k.inv(lc);
            a.iter().map(|c| k.mul(c, &inv)).collect()
        }
    }
}

/// `(quotient, remainder)`.
pub fn divmod(k: &Ext, a: &Poly, b: &Poly) -> (Poly, Poly) {
    assert!(!b.is_empty(), "division by the zero polynomial");
    let db = deg(b) as usize;
    let binv = k.inv(b.last().unwrap());
    let mut r = a.clone();
    let mut q = vec![k.zero(); a.len().saturating_sub(db)];
    loop {
        r = norm(k, r);
        if deg(&r) < db as isize {
            break;
        }
        let dr = deg(&r) as usize;
        let c = k.mul(r.last().unwrap(), &binv);
        q[dr - db] = c.clone();
        for (j, bj) in b.iter().enumerate() {
            let t = k.mul(&c, bj);
            r[dr - db + j] = k.sub(&r[dr - db + j], &t);
        }
    }
    (norm(k, q), r)
}

pub fn rem(k: &Ext, a: &Poly, b: &Poly) -> Poly {
    divmod(k, a, b).1
}

pub fn gcd(k: &Ext, a: &Poly, b: &Poly) -> Poly {
    let (mut a, mut b) = (norm(k, a.clone()), norm(k, b.clone()));
    while !b.is_empty() {
        let t = rem(k, &a, &b);
        a = b;
        b = t;
    }
    monic(k, &a)
}

/// `base^e mod f`, with `e` a little-endian limb vector.
pub fn powmod(k: &Ext, base: &Poly, e: &[u64], f: &Poly) -> Poly {
    let bl = bn_bitlen(e);
    let mut r: Poly = vec![k.one()];
    let b = rem(k, base, f);
    for i in (0..bl).rev() {
        r = rem(k, &mul(k, &r, &r), f);
        if bn_bit(e, i) {
            r = rem(k, &mul(k, &r, &b), f);
        }
    }
    r
}

pub fn evaluate(k: &Ext, f: &Poly, x: &Elem) -> Elem {
    let mut acc = k.zero();
    for c in f.iter().rev() {
        acc = k.mul(&acc, x);
        acc = k.add(&acc, c);
    }
    acc
}

// ------------------------------------------------------- Frobenius ladder

/// Given `A = X^(q^a) mod F` and `B = X^(q^b) mod F`, return `X^(q^(a+b))`.
///
/// `B(X)^(q^a) = sum_i B_i^(q^a) * (X^(q^a))^i`, and `X^(q^a) = A mod F`.
fn compose_frob(k: &Ext, f: &Poly, a: &Poly, ea: usize, b: &Poly) -> Poly {
    if b.is_empty() {
        return Vec::new();
    }
    let mut pw: Vec<Poly> = vec![vec![k.one()]];
    for _ in 0..deg(b) {
        let next = rem(k, &mul(k, pw.last().unwrap(), a), f);
        pw.push(next);
    }
    let mut out: Poly = Vec::new();
    for (i, c) in b.iter().enumerate() {
        if !k.is_zero(c) {
            let s = k.frob(c, ea);
            out = add(k, &out, &scal(k, &pw[i], &s));
        }
    }
    out
}

/// `X^(q^n) mod F` by iterating the Frobenius `n-1` times.
///
/// Costs `O(n * D^2)` products against the ladder's `O(log n * D^3)`, so it
/// wins whenever `D` is large relative to `n` -- which is exactly the `q = 13`
/// and `q = 23` parameter sets, where `D` reaches 24. [`x_pow_qn`] picks.
fn x_pow_qn_linear(k: &Ext, f: &Poly) -> Poly {
    let d = deg(f).max(0) as usize;
    let mut xq: Poly = vec![k.zero(); k.q as usize];
    xq.push(k.one());
    let q = rem(k, &xq, f);
    // Q^i mod f, computed once.
    let mut qpow: Vec<Poly> = vec![vec![k.one()]];
    for _ in 0..d {
        let next = rem(k, &mul(k, qpow.last().unwrap(), &q), f);
        qpow.push(next);
    }
    let mut cur = q;
    for _ in 1..k.n {
        // A_{j+1} = sum_i (A_j)_i^q * Q^i, since the coefficients are the
        // only part the Frobenius moves.
        let mut out: Poly = Vec::new();
        for (i, c) in cur.iter().enumerate() {
            if !k.is_zero(c) {
                let s = k.frob(c, 1);
                out = add(k, &out, &scal(k, &qpow[i], &s));
            }
        }
        cur = out;
    }
    cur
}

/// `X^(q^n) mod F`, by whichever Frobenius schedule is cheaper here.
pub fn x_pow_qn(k: &Ext, f: &Poly) -> Poly {
    let d = deg(f).max(0) as usize;
    let bits = (usize::BITS - k.n.leading_zeros()) as usize;
    // linear ~ n*D^2 products, ladder ~ 2*bits*D^3
    if k.n < 2 * bits * d {
        return x_pow_qn_linear(k, f);
    }
    x_pow_qn_ladder(k, f)
}

/// `X^(q^n) mod F`, by doubling the Frobenius exponent rather than the power.
fn x_pow_qn_ladder(k: &Ext, f: &Poly) -> Poly {
    let n = k.n;
    let mut xq: Poly = vec![k.zero(); k.q as usize];
    xq.push(k.one()); // X^q
    let q = rem(k, &xq, f);
    let mut res: Option<Poly> = None;
    let mut exp = 0usize;
    for i in (0..usize::BITS - n.leading_zeros()).rev() {
        let bit = (n >> i) & 1;
        match res {
            None => {
                res = Some(q.clone());
                exp = 1;
            }
            Some(cur) => {
                let mut t = compose_frob(k, f, &cur, exp, &cur);
                exp *= 2;
                if bit == 1 {
                    t = compose_frob(k, f, &q, 1, &t);
                    exp += 1;
                }
                res = Some(t);
            }
        }
    }
    debug_assert_eq!(exp, n);
    res.unwrap_or_default()
}

// ------------------------------------------------------------ root finding

fn rand_nonzero(k: &Ext, xof: &mut Xof) -> Elem {
    loop {
        let cs = sample_fq(xof, k.n, k.q);
        let e = k.from_coords(&cs);
        if !k.is_zero(&e) {
            return e;
        }
    }
}

/// Split a product of distinct linear factors into its roots.
fn split(k: &Ext, g: &Poly, xof: &mut Xof, out: &mut Vec<Elem>) {
    if deg(g) == 1 {
        out.push(k.neg(&g[0])); // g is monic: X + g0
        return;
    }
    loop {
        let delta = rand_nonzero(k, xof);
        let c = if k.fq.p == 2 {
            // characteristic 2: the absolute trace to F_2 separates the roots
            let m = k.n * k.fq.k as usize;
            let mut h: Poly = Vec::new();
            let mut term: Poly = vec![k.zero(), delta];
            for _ in 0..m {
                h = add(k, &h, &term);
                term = rem(k, &mul(k, &term, &term), g);
            }
            gcd(k, g, &h)
        } else {
            let mut e = bn_pow(k.q, k.n);
            bn_sub1(&mut e);
            bn_divmod_small(&mut e, 2); // (q^n - 1) / 2
            let base: Poly = vec![delta, k.one()];
            let h = powmod(k, &base, &e, g);
            gcd(k, g, &sub(k, &h, &vec![k.one()]))
        };
        if deg(&c) > 0 && deg(&c) < deg(g) {
            let (other, r) = divmod(k, g, &c);
            debug_assert!(r.is_empty());
            split(k, &monic(k, &c), xof, out);
            split(k, &monic(k, &other), xof, out);
            return;
        }
    }
}

/// Every root of `f` in `K`, in canonical order.
///
/// `gcd(X^(q^n) - X, f)` strips `f` to the product of its distinct linear
/// factors; equal-degree splitting peels those apart.
///
/// The result is sorted by coordinate vector read as a base-`q` integer with
/// `c_{n-1}` most significant. That ordering is normative: splitting returns
/// roots in an order that depends on the random `delta` it happened to draw,
/// and the signer takes the first root passing the IP check, so without a
/// canonical order two conforming implementations could emit different -- both
/// valid -- signatures for the same key and message.
pub fn roots(k: &Ext, f: &Poly, xof: &mut Xof) -> Vec<Elem> {
    let f = monic(k, &norm(k, f.clone()));
    if deg(&f) < 1 {
        return Vec::new();
    }
    if deg(&f) == 1 {
        return vec![k.neg(&f[0])];
    }
    let a = x_pow_qn(k, &f);
    let x: Poly = vec![k.zero(), k.one()];
    let g = monic(k, &gcd(k, &f, &sub(k, &a, &x)));
    if deg(&g) < 1 {
        return Vec::new();
    }
    let mut out = Vec::new();
    split(k, &g, xof, &mut out);
    out.sort_by(|p, s| {
        let (mut cp, mut cs) = (k.coords(p), k.coords(s));
        cp.reverse();
        cs.reverse();
        cp.cmp(&cs)
    });
    out
}

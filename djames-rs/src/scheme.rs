//! Key generation, signing and verification (`d-james-spec.md` §5-§8).
//!
//! Row-vector convention: the signature is `a`, the hidden HFE variable is
//! `x = a S`, and `a = x S^-1`.

use crate::codec::{self, CodecError, DigitPacker};
use crate::ct;
use crate::gf::{Elem, Ext, Fq, PackedVec};
use crate::linalg::{self, Matrix};
use crate::params::Params;
use crate::poly::{self, Poly};
use crate::symmetric::{hash_to_fq, sample_fq, Xof, DOM_EDF, DOM_KEY, DOM_PAD, DOM_PRF};
use alloc::vec;
use alloc::vec::Vec;

/// A secret key is exactly its seed, so its length is part of the interface.
pub const SEED_BYTES: usize = 32;

/// The signer walks `salt = 0, 1, 2, ...` and the verifier re-walks it,
/// because the salt is not transmitted.
///
/// Security-relevant, not merely an engineering limit: a forger may aim at any
/// of these hash values, conceding `log2(MAX_SALT) = 8` bits. With the paper's
/// expected trial counts (1.582 in characteristic 2, 3.07 otherwise) 256 salts
/// fail with probability about `2^-369` resp. `2^-146`.
pub const MAX_SALT: u64 = 256;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Error {
    /// The seed was not [`SEED_BYTES`] long.
    SeedLength,
    /// Signing exhausted [`MAX_SALT`]; astronomically unlikely.
    SigningFailed,
    /// A serialized value was malformed.
    Codec(CodecError),
}

impl From<CodecError> for Error {
    fn from(e: CodecError) -> Self {
        Error::Codec(e)
    }
}

pub struct SecretKey {
    pub params: &'static Params,
    seed: [u8; SEED_BYTES],
    prf: [u8; 32],
    k: Ext,
    lam: Vec<Elem>,
    mbil: Vec<Vec<Elem>>,
    gq: Vec<Vec<Elem>>,
    dragon: Vec<Vec<Elem>>,
    mz: Matrix,
    s: Matrix,
    sinv: Matrix,
    tinv: Matrix,
}

impl SecretKey {
    /// The seed, which is the whole secret key.
    pub fn to_bytes(&self) -> [u8; SEED_BYTES] {
        self.seed
    }
}

impl Drop for SecretKey {
    fn drop(&mut self) {
        ct::zeroize8(&mut self.seed);
        ct::zeroize8(&mut self.prf);
        for e in self.lam.iter_mut() {
            ct::zeroize64(e);
        }
        for row in self.s.iter_mut().chain(self.sinv.iter_mut()) {
            ct::zeroize8(row);
        }
    }
}

/// `m` quadratic equations over `F_q`, stored transposed: one packed `F_q^m`
/// vector per monomial, holding that monomial's coefficient in all `m`
/// equations at once. Evaluation is then a handful of packed additions.
pub struct PublicKey {
    pub params: &'static Params,
    w: PackedVec,
    aa: Vec<Vec<u64>>,
    ab: Vec<Vec<u64>>,
    off: Vec<usize>,
}

fn offsets(n: usize) -> Vec<usize> {
    let mut off = Vec::with_capacity(n);
    let mut acc = 0;
    for i in 0..n {
        off.push(acc);
        acc += n - i;
    }
    off
}

// ------------------------------------------------------------- key  gen

pub fn keygen(p: &'static Params, seed: &[u8]) -> Result<(PublicKey, SecretKey), Error> {
    if seed.len() != SEED_BYTES {
        return Err(Error::SeedLength);
    }
    let tag = p.tag();
    let mut xof = Xof::new(&[DOM_KEY, seed, tag.as_bytes()]);
    let k = Ext::new(p.q, p.n, p.fpoly);
    let (n, r, d, ny) = (p.n, p.r, p.d, p.ny);

    // Draw order is normative; see the spec.
    let mut lam = Vec::with_capacity(p.monomials.len());
    for _ in 0..p.monomials.len() {
        lam.push(sample_k_nonzero(&k, &mut xof));
    }
    let mut mbil = Vec::with_capacity(d);
    for _ in 0..d {
        mbil.push((0..r).map(|_| sample_k(&k, &mut xof)).collect::<Vec<_>>());
    }
    let mut gq: Vec<Vec<Elem>> = (0..r).map(|_| vec![k.zero(); r]).collect();
    for i in 0..r {
        for j in i..r {
            gq[i][j] = sample_k(&k, &mut xof);
        }
    }
    let mut dragon: Vec<Vec<Elem>> = Vec::new();
    if p.is_dragon() {
        for _ in 0..=d {
            dragon.push((0..ny).map(|_| sample_k(&k, &mut xof)).collect());
        }
    }
    let mz = linalg::random_full_rank(&mut xof, r, n, p.q);
    let (s, sinv) = linalg::random_invertible(&mut xof, n, p.q);
    let (t, tinv) = linalg::random_invertible(&mut xof, n, p.q);

    let mut prf = [0u8; 32];
    Xof::new(&[DOM_PRF, seed, tag.as_bytes()]).read(&mut prf);

    let pk = build_public_key(p, &k, &lam, &mbil, &gq, &dragon, &mz, &s, &t);

    let mut seed_arr = [0u8; SEED_BYTES];
    seed_arr.copy_from_slice(seed);
    let sk = SecretKey {
        params: p,
        seed: seed_arr,
        prf,
        k,
        lam,
        mbil,
        gq,
        dragon,
        mz,
        s,
        sinv,
        tinv,
    };
    Ok((pk, sk))
}

fn sample_k(k: &Ext, xof: &mut Xof) -> Elem {
    let cs = sample_fq(xof, k.n, k.q);
    k.from_coords(&cs)
}

fn sample_k_nonzero(k: &Ext, xof: &mut Xof) -> Elem {
    loop {
        let e = sample_k(k, xof);
        if !k.is_zero(&e) {
            return e;
        }
    }
}

/// One rank-one outer product `u (x) v` over K.
type Term = (Vec<Elem>, Vec<Elem>);

#[allow(clippy::too_many_arguments)]
fn build_public_key(
    p: &'static Params,
    k: &Ext,
    lam: &[Elem],
    mbil: &[Vec<Elem>],
    gq: &[Vec<Elem>],
    dragon: &[Vec<Elem>],
    mz: &Matrix,
    s: &Matrix,
    t: &Matrix,
) -> PublicKey {
    let fq = Fq::new(p.q);
    let (n, m, r, d, ny) = (p.n, p.m, p.r, p.d, p.ny);

    // The central map's inputs, as linear forms in (a, b) over K.
    let sigma: Vec<Elem> = (0..n).map(|i| k.from_coords(&s[i])).collect();
    let a_pow: Vec<Vec<Elem>> = (0..=d)
        .map(|kk| sigma.iter().map(|x| k.frob(x, kk)).collect())
        .collect();
    let mut zf: Vec<Vec<u8>> = Vec::with_capacity(r);
    for j in 0..r {
        let mut col = Vec::with_capacity(n);
        for i in 0..n {
            let mut acc = 0u8;
            for c in 0..n {
                acc = fq.add(acc, fq.mul(mz[j][c], s[i][c]));
            }
            col.push(acc);
        }
        zf.push(col);
    }
    let one = k.one();
    let zk: Vec<Vec<Elem>> = zf
        .iter()
        .map(|col| col.iter().map(|&c| k.scal(&one, c)).collect())
        .collect();

    // Every piece of H is a product of two linear forms, hence a rank-one
    // outer product. Folding the IP families over their inner index first
    // cuts the term count from d*r + r(r+1)/2 to d + r.
    let mut aa_terms: Vec<Term> = Vec::new();
    for (ti, &(mi, mj)) in p.monomials.iter().enumerate() {
        let u: Vec<Elem> = a_pow[mi].iter().map(|x| k.mul(&lam[ti], x)).collect();
        aa_terms.push((u, a_pow[mj].clone()));
    }
    for kk in 0..d {
        let mut v = vec![k.zero(); n];
        for j in 0..r {
            if !k.is_zero(&mbil[kk][j]) {
                for i in 0..n {
                    let t = k.scal(&mbil[kk][j], zf[j][i]);
                    v[i] = k.add(&v[i], &t);
                }
            }
        }
        aa_terms.push((a_pow[kk].clone(), v));
    }
    for i0 in 0..r {
        let mut v = vec![k.zero(); n];
        for j0 in i0..r {
            if !k.is_zero(&gq[i0][j0]) {
                for i in 0..n {
                    let t = k.scal(&gq[i0][j0], zf[j0][i]);
                    v[i] = k.add(&v[i], &t);
                }
            }
        }
        aa_terms.push((zk[i0].clone(), v));
    }
    let mut ab_terms: Vec<Term> = Vec::new();
    if p.is_dragon() {
        for kk in 0..=d {
            ab_terms.push((a_pow[kk].clone(), dragon[kk].clone()));
        }
    }

    let off = offsets(n);
    let mut sc = k.scratch();
    let mut prod = k.zero();
    let mut c_aa: Vec<Elem> = vec![k.zero(); n * (n + 1) / 2];
    for (u, v) in &aa_terms {
        for pp in 0..n {
            if k.is_zero(&u[pp]) {
                continue;
            }
            for ss in 0..n {
                k.mul_into(&mut prod, &u[pp], &v[ss], &mut sc);
                let idx = if pp <= ss {
                    off[pp] + ss - pp
                } else {
                    off[ss] + pp - ss
                };
                k.add_assign(&mut c_aa[idx], &prod);
            }
        }
    }
    let mut c_ab: Vec<Elem> = if ny > 0 {
        vec![k.zero(); n * ny]
    } else {
        Vec::new()
    };
    for (u, v) in &ab_terms {
        for pp in 0..n {
            if k.is_zero(&u[pp]) {
                continue;
            }
            let base = pp * ny;
            for ss in 0..ny {
                k.mul_into(&mut prod, &u[pp], &v[ss], &mut sc);
                k.add_assign(&mut c_ab[base + ss], &prod);
            }
        }
    }

    // Project the n central equations onto the m published ones: the minus
    // modifier is exactly "keep the first m rows of T".
    let w = PackedVec::new(p.q, m);
    let pi: Vec<Vec<u64>> = (0..n)
        .map(|c| {
            let col: Vec<u8> = (0..m).map(|row| t[row][c]).collect();
            w.from_digits(&col)
        })
        .collect();
    let project = |e: &Elem| -> Vec<u64> {
        let mut acc = w.zero();
        k.project_into(&mut acc, e, &pi, &w);
        acc
    };
    let aa: Vec<Vec<u64>> = c_aa.iter().map(&project).collect();
    let ab: Vec<Vec<u64>> = c_ab.iter().map(&project).collect();
    PublicKey {
        params: p,
        w,
        aa,
        ab,
        off,
    }
}

impl PublicKey {
    /// The `m` public equations at `(a, b)`, as one packed `F_q^m` vector.
    ///
    /// Everything here is public, so the zero-skipping is free of timing
    /// concerns and worth the speed.
    fn eval(&self, a: &[u8], b: Option<&[u8]>) -> Vec<u64> {
        let p = self.params;
        let fq = Fq::new(p.q);
        let mut acc = self.w.zero();
        for pp in 0..p.n {
            if a[pp] == 0 {
                continue;
            }
            for ss in pp..p.n {
                let c = fq.mul(a[pp], a[ss]);
                if c != 0 {
                    self.w
                        .scal_add_assign(&mut acc, &self.aa[self.off[pp] + ss - pp], c);
                }
            }
        }
        if let Some(b) = b {
            for pp in 0..p.n {
                if a[pp] == 0 {
                    continue;
                }
                for ss in 0..p.ny {
                    let c = fq.mul(a[pp], b[ss]);
                    if c != 0 {
                        self.w
                            .scal_add_assign(&mut acc, &self.ab[pp * p.ny + ss], c);
                    }
                }
            }
        }
        acc
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut pk = DigitPacker::new(self.params.q);
        for v in self.aa.iter().chain(self.ab.iter()) {
            pk.push(&self.w.to_digits(v));
        }
        pk.finish()
    }

    pub fn from_bytes(p: &'static Params, data: &[u8]) -> Result<Self, Error> {
        let w = PackedVec::new(p.q, p.m);
        let total = p.pk_coeffs() * p.m;
        let digits = codec::bytes_to_digits(data, total, p.q)?;
        let mut vecs: Vec<Vec<u64>> = digits.chunks(p.m).map(|c| w.from_digits(c)).collect();
        let n_aa = p.n * (p.n + 1) / 2;
        let ab = vecs.split_off(n_aa);
        Ok(PublicKey {
            params: p,
            w,
            aa: vecs,
            ab,
            off: offsets(p.n),
        })
    }
}

// ---------------------------------------------------------------- signing

/// The univariate polynomial to root-find, for one hash and one IP guess.
fn central_poly(p: &Params, sk: &SecretKey, consts: Option<&[Elem]>, z: &[u8]) -> Poly {
    let k = &sk.k;
    let fq = Fq::new(p.q);
    let q = p.q as usize;
    let mut f: Poly = vec![k.zero(); p.dd + 1];
    for (ti, &(i, j)) in p.monomials.iter().enumerate() {
        let e = q.pow(i as u32) + q.pow(j as u32);
        f[e] = k.add(&f[e], &sk.lam[ti]);
    }
    if let Some(c) = consts {
        for kk in 0..=p.d {
            let e = q.pow(kk as u32);
            f[e] = k.add(&f[e], &c[kk]);
        }
    }
    for kk in 0..p.d {
        let mut acc = k.zero();
        for j in 0..p.r {
            if z[j] != 0 {
                let t = k.scal(&sk.mbil[kk][j], z[j]);
                acc = k.add(&acc, &t);
            }
        }
        let e = q.pow(kk as u32);
        f[e] = k.add(&f[e], &acc);
    }
    let mut cst = k.zero();
    for i in 0..p.r {
        for j in i..p.r {
            let c = fq.mul(z[i], z[j]);
            if c != 0 {
                let t = k.scal(&sk.gq[i][j], c);
                cst = k.add(&cst, &t);
            }
        }
    }
    f[0] = k.add(&f[0], &cst);
    f
}

/// Sign `msg`. Deterministic: no entropy beyond the key and the message.
pub fn sign(p: &'static Params, sk: &SecretKey, msg: &[u8]) -> Result<Vec<u8>, Error> {
    let k = &sk.k;
    let fq = Fq::new(p.q);
    let (r, q) = (p.r, p.q);
    let mut edf = Xof::new(&[DOM_EDF, &sk.prf, msg]);
    let nz = (q as u64).pow(r as u32);

    for counter in 0..MAX_SALT {
        let mut consts: Vec<Elem> = Vec::new();
        let mut target: Option<Elem> = None;
        if p.is_dragon() {
            let y = hash_to_fq(msg, counter, p.ny, q);
            for kk in 0..=p.d {
                let mut acc = k.zero();
                for j in 0..p.ny {
                    if y[j] != 0 {
                        let t = k.scal(&sk.dragon[kk][j], y[j]);
                        acc = k.add(&acc, &t);
                    }
                }
                consts.push(acc);
            }
        } else {
            let h = hash_to_fq(msg, 0, p.m, q);
            let mut px = Xof::new(&[DOM_PAD, &sk.prf, msg, &counter.to_le_bytes()]);
            let pad = sample_fq(&mut px, p.a, q);
            let mut c = h;
            c.extend_from_slice(&pad);
            let u = linalg::mat_vec(&fq, &sk.tinv, &c);
            target = Some(k.from_coords(&u));
        }

        for e in 0..nz {
            let z: Vec<u8> = (0..r)
                .map(|j| ((e / (q as u64).pow((r - 1 - j) as u32)) % q as u64) as u8)
                .collect();
            let mut f = central_poly(p, sk, if p.is_dragon() { Some(&consts) } else { None }, &z);
            if let Some(t) = &target {
                f[0] = k.sub(&f[0], t);
            }
            let f = poly::norm(k, f);
            for root in poly::roots(k, &f, &mut edf) {
                if k.is_zero(&root) {
                    continue; // homogeneous system: 0 is the trivial root
                }
                let x = k.coords(&root);
                if linalg::mat_vec(&fq, &sk.mz, &x) != z {
                    continue; // the IP guess was wrong
                }
                let a = linalg::vec_mat(&fq, &x, &sk.sinv);
                if a.iter().any(|&c| c != 0) {
                    return Ok(codec::encode_vector(&a, q));
                }
            }
        }
    }
    Err(Error::SigningFailed)
}

// ----------------------------------------------------------- verification

pub fn verify(p: &'static Params, pk: &PublicKey, msg: &[u8], sig: &[u8]) -> bool {
    let a = match codec::decode_vector(sig, p.n, p.q) {
        Ok(a) => a,
        Err(_) => return false,
    };
    if a.iter().all(|&c| c == 0) {
        return false; // the trivial root of a homogeneous system
    }
    if p.is_dragon() {
        for salt in 0..MAX_SALT {
            let b = hash_to_fq(msg, salt, p.ny, p.q);
            if pk.w.is_zero(&pk.eval(&a, Some(&b))) {
                return true;
            }
        }
        false
    } else {
        let h = hash_to_fq(msg, 0, p.m, p.q);
        pk.eval(&a, None) == pk.w.from_digits(&h)
    }
}

/// SHA3-256 of the serialized public key, as used by the test vectors.
pub fn pk_digest(pk: &PublicKey) -> [u8; 32] {
    crate::keccak::sha3_256(&pk.to_bytes())
}

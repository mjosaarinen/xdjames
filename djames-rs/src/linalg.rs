//! Matrices over `F_q`: sampling, inversion, rank.
//!
//! Rows are `Vec<u8>`, one coefficient per byte, row-major.
//!
//! # Constant time
//!
//! Elimination runs a fixed nest of loops over every row and column and uses
//! masks rather than branches, so for a **full-rank square** input -- which is
//! what `S` and `T` always are -- the path taken is independent of the
//! entries. The rank counter advances on a mask, so a rank-deficient input
//! does take a different path; that only ever happens to a *candidate* matrix
//! that is then rejected and discarded, never to accepted key material.
//!
//! The one residual leak is `random_full_rank` for `M_Z` (`r x n`, `r = 2`):
//! which columns end up as pivots depends on the accepted matrix. For a
//! uniform matrix these are the first `r` columns with overwhelming
//! probability, so the leak is slight, but it is real and worth naming.

use crate::ct;
use crate::gf::Fq;
use crate::symmetric::{sample_fq, Xof};
use alloc::vec;
use alloc::vec::Vec;

pub type Matrix = Vec<Vec<u8>>;

/// A uniform `rows x cols` matrix, filled row-major from the stream.
pub fn random_matrix(xof: &mut Xof, rows: usize, cols: usize, q: u32) -> Matrix {
    let flat = sample_fq(xof, rows * cols, q);
    flat.chunks(cols).map(|c| c.to_vec()).collect()
}

pub fn identity(n: usize) -> Matrix {
    (0..n)
        .map(|i| (0..n).map(|j| u8::from(i == j)).collect())
        .collect()
}

/// Row-reduce `m` (carrying `aug` along), returning the rank.
fn echelon(fq: &Fq, m: &mut Matrix, ncols: usize, mut aug: Option<&mut Matrix>) -> usize {
    let nrows = m.len();
    let mut rank = 0usize;
    for c in 0..ncols {
        if rank >= nrows {
            break;
        }
        // Fold every lower row into the pivot row under the mask "the pivot
        // row's entry is still zero". After this the pivot is nonzero iff the
        // column has any nonzero entry at or below `rank`.
        for i in rank + 1..nrows {
            let need = ct::eq_mask32(m[rank][c] as u32, 0) as u8 & 1;
            for j in 0..ncols {
                let v = fq.mul(m[i][j], need);
                m[rank][j] = fq.add(m[rank][j], v);
            }
            if let Some(a) = aug.as_deref_mut() {
                for j in 0..a[i].len() {
                    let v = fq.mul(a[i][j], need);
                    a[rank][j] = fq.add(a[rank][j], v);
                }
            }
        }
        let piv = m[rank][c];
        if piv == 0 {
            continue; // rank-deficient: only reachable for rejected candidates
        }
        let inv = fq.inv(piv);
        for j in 0..ncols {
            m[rank][j] = fq.mul(m[rank][j], inv);
        }
        if let Some(a) = aug.as_deref_mut() {
            for j in 0..a[rank].len() {
                a[rank][j] = fq.mul(a[rank][j], inv);
            }
        }
        for i in 0..nrows {
            if i == rank {
                continue;
            }
            let f = fq.neg(m[i][c]);
            if f != 0 {
                for j in 0..ncols {
                    let v = fq.mul(m[rank][j], f);
                    m[i][j] = fq.add(m[i][j], v);
                }
                if let Some(a) = aug.as_deref_mut() {
                    for j in 0..a[rank].len() {
                        let v = fq.mul(a[rank][j], f);
                        a[i][j] = fq.add(a[i][j], v);
                    }
                }
            }
        }
        rank += 1;
    }
    rank
}

pub fn rank(fq: &Fq, m: &Matrix, ncols: usize) -> usize {
    let mut work = m.clone();
    echelon(fq, &mut work, ncols, None)
}

/// Inverse of a square matrix, or `None` if singular.
pub fn invert(fq: &Fq, m: &Matrix) -> Option<Matrix> {
    let n = m.len();
    let mut work = m.clone();
    let mut aug = identity(n);
    let rk = echelon(fq, &mut work, n, Some(&mut aug));
    if rk == n {
        Some(aug)
    } else {
        None
    }
}

/// A uniform invertible `n x n` matrix and its inverse.
///
/// Rejection is over the *sampled* matrix, so the retry count is independent
/// of the value finally returned. A uniform matrix over `F_q` is invertible
/// with probability `prod (1 - q^-i)`, about 0.29 for `q = 2` and rising with
/// `q`, so the loop turns over a handful of times at worst.
pub fn random_invertible(xof: &mut Xof, n: usize, q: u32) -> (Matrix, Matrix) {
    let fq = Fq::new(q);
    loop {
        let m = random_matrix(xof, n, n, q);
        if let Some(inv) = invert(&fq, &m) {
            return (m, inv);
        }
    }
}

/// A uniform `rows x cols` matrix of full row rank.
pub fn random_full_rank(xof: &mut Xof, rows: usize, cols: usize, q: u32) -> Matrix {
    let fq = Fq::new(q);
    loop {
        let m = random_matrix(xof, rows, cols, q);
        if rank(&fq, &m, cols) == rows {
            return m;
        }
    }
}

/// `M . v`, with `v` a column of length `cols`.
pub fn mat_vec(fq: &Fq, m: &Matrix, v: &[u8]) -> Vec<u8> {
    m.iter()
        .map(|row| {
            let mut acc = 0u8;
            for (&c, &x) in row.iter().zip(v.iter()) {
                acc = fq.add(acc, fq.mul(c, x));
            }
            acc
        })
        .collect()
}

/// `v . M`, with `v` a row of length `rows`.
pub fn vec_mat(fq: &Fq, v: &[u8], m: &Matrix) -> Vec<u8> {
    let cols = m[0].len();
    let mut out = vec![0u8; cols];
    for (&x, row) in v.iter().zip(m.iter()) {
        for j in 0..cols {
            out[j] = fq.add(out[j], fq.mul(x, row[j]));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symmetric::Xof;

    #[test]
    fn inverse_round_trips() {
        for q in [2u32, 4, 5, 13, 23] {
            let fq = Fq::new(q);
            let mut x = Xof::new(&[b"la", &[q as u8]]);
            let (m, mi) = random_invertible(&mut x, 12, q);
            for i in 0..12 {
                let e: Vec<u8> = (0..12).map(|k| u8::from(k == i)).collect();
                assert_eq!(vec_mat(&fq, &vec_mat(&fq, &e, &m), &mi), e, "q={q}");
                assert_eq!(vec_mat(&fq, &vec_mat(&fq, &e, &mi), &m), e, "q={q}");
            }
            let mut x = Xof::new(&[b"z", &[q as u8]]);
            let z = random_full_rank(&mut x, 3, 12, q);
            assert_eq!(rank(&fq, &z, 12), 3);
        }
    }
}

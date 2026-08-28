//! Parameter sets and canonical field polynomials.
//!
//! Generated from `djames-py` by `tools/gen_params_rs.py`; do not edit by
//! hand. The tables mirror `d-james-spec.md` §2 and §3.2 exactly, so the two
//! implementations cannot drift apart through a transcription slip.

use alloc::string::String;
use alloc::vec::Vec;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Scheme {
    /// HFE-_IP: verification checks `P(a) = Hash(msg)`.
    James,
    /// HFE-_IP with Dragon terms: the system is homogeneous, `P(a, b) = 0`.
    DJames,
}

/// One parameter set. `m = n - a` and `D = q^d + 1` always.
#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub name: &'static str,
    pub scheme: Scheme,
    pub q: u32,
    pub n: usize,
    pub a: usize,
    pub m: usize,
    /// Hash length in F_q symbols; 0 for James.
    pub ny: usize,
    pub r: usize,
    pub d: usize,
    pub dd: usize,
    /// `(i, j)` with a nonzero `lambda_{i,j} X^(q^i + q^j)`.
    pub monomials: &'static [(usize, usize)],
    /// `f = t^n + sum c_i t^i`, low coefficients.
    pub fpoly: &'static [u8],
}

impl Params {
    pub fn is_dragon(&self) -> bool {
        matches!(self.scheme, Scheme::DJames)
    }

    /// Number of published quadratic coefficients.
    pub fn pk_coeffs(&self) -> usize {
        self.n * (self.n + 1) / 2 + self.n * self.ny
    }

    /// Exact serialized signature length.
    pub fn sig_bytes(&self) -> usize {
        crate::codec::vec_bytes(self.n, self.q)
    }

    /// The domain-separation tag of `d-james-spec.md` §2.3.
    pub fn tag(&self) -> String {
        use core::fmt::Write;
        let mut s = String::new();
        let scheme = if self.is_dragon() { "d-james" } else { "james" };
        write!(
            s,
            "D-James/v1/{}/q{}/n{}/a{}/r{}/D{}/ny{}/mon",
            scheme, self.q, self.n, self.a, self.r, self.dd, self.ny
        )
        .unwrap();
        for (k, (i, j)) in self.monomials.iter().enumerate() {
            if k > 0 {
                s.push('.');
            }
            write!(s, "{i}-{j}").unwrap();
        }
        s
    }
}

#[rustfmt::skip]
static FP_2_32: &[u8] = &[1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_2_189: &[u8] = &[1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_2_283: &[u8] = &[1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_2_390: &[u8] = &[1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_2_578: &[u8] = &[1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_4_24: &[u8] = &[2, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_4_105: &[u8] = &[1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_4_149: &[u8] = &[2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_4_223: &[u8] = &[3, 0, 3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_4_309: &[u8] = &[2, 0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_5_21: &[u8] = &[1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_5_94: &[u8] = &[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_5_132: &[u8] = &[1, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_5_206: &[u8] = &[3, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_5_274: &[u8] = &[1, 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_13_18: &[u8] = &[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_13_77: &[u8] = &[2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_13_91: &[u8] = &[8, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_13_171: &[u8] = &[3, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_13_192: &[u8] = &[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_23_16: &[u8] = &[20, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_23_74: &[u8] = &[3, 9, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_23_78: &[u8] = &[3, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_23_163: &[u8] = &[3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
#[rustfmt::skip]
static FP_23_167: &[u8] = &[7, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

static MON_D_JAMES_128_Q13: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_D_JAMES_128_Q2: &[(usize, usize)] = &[(0, 1), (0, 2)];
static MON_D_JAMES_128_Q23: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_D_JAMES_128_Q4: &[(usize, usize)] = &[(0, 0), (0, 2)];
static MON_D_JAMES_128_Q5: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_D_JAMES_256_Q13: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_D_JAMES_256_Q2: &[(usize, usize)] = &[(0, 1), (0, 2)];
static MON_D_JAMES_256_Q23: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_D_JAMES_256_Q4: &[(usize, usize)] = &[(0, 0), (0, 2)];
static MON_D_JAMES_256_Q5: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_JAMES_128_Q13: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_JAMES_128_Q2: &[(usize, usize)] = &[(0, 1), (0, 2)];
static MON_JAMES_128_Q23: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_JAMES_128_Q4: &[(usize, usize)] = &[(0, 0), (0, 2)];
static MON_JAMES_128_Q5: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_JAMES_256_Q13: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_JAMES_256_Q2: &[(usize, usize)] = &[(0, 1), (0, 2)];
static MON_JAMES_256_Q23: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_JAMES_256_Q4: &[(usize, usize)] = &[(0, 0), (0, 2)];
static MON_JAMES_256_Q5: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_TOY_D_JAMES_Q13: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_TOY_D_JAMES_Q2: &[(usize, usize)] = &[(0, 1), (0, 2)];
static MON_TOY_D_JAMES_Q23: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_TOY_D_JAMES_Q4: &[(usize, usize)] = &[(0, 0), (0, 2)];
static MON_TOY_D_JAMES_Q5: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_TOY_JAMES_Q13: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_TOY_JAMES_Q2: &[(usize, usize)] = &[(0, 1), (0, 2)];
static MON_TOY_JAMES_Q23: &[(usize, usize)] = &[(0, 0), (0, 1)];
static MON_TOY_JAMES_Q4: &[(usize, usize)] = &[(0, 0), (0, 2)];
static MON_TOY_JAMES_Q5: &[(usize, usize)] = &[(0, 0), (0, 1)];

/// Every parameter set, in the order of `d-james-spec.md` §2.
pub static ALL: &[Params] = &[
    Params {
        name: "d-james-128-q13",
        scheme: Scheme::DJames,
        q: 13,
        n: 77,
        a: 21,
        m: 56,
        ny: 70,
        r: 2,
        d: 1,
        dd: 14,
        monomials: MON_D_JAMES_128_Q13,
        fpoly: FP_13_77,
    },
    Params {
        name: "d-james-128-q2",
        scheme: Scheme::DJames,
        q: 2,
        n: 189,
        a: 27,
        m: 162,
        ny: 256,
        r: 2,
        d: 2,
        dd: 5,
        monomials: MON_D_JAMES_128_Q2,
        fpoly: FP_2_189,
    },
    Params {
        name: "d-james-128-q23",
        scheme: Scheme::DJames,
        q: 23,
        n: 74,
        a: 21,
        m: 53,
        ny: 57,
        r: 2,
        d: 1,
        dd: 24,
        monomials: MON_D_JAMES_128_Q23,
        fpoly: FP_23_74,
    },
    Params {
        name: "d-james-128-q4",
        scheme: Scheme::DJames,
        q: 4,
        n: 105,
        a: 21,
        m: 84,
        ny: 128,
        r: 2,
        d: 2,
        dd: 17,
        monomials: MON_D_JAMES_128_Q4,
        fpoly: FP_4_105,
    },
    Params {
        name: "d-james-128-q5",
        scheme: Scheme::DJames,
        q: 5,
        n: 94,
        a: 21,
        m: 73,
        ny: 111,
        r: 2,
        d: 1,
        dd: 6,
        monomials: MON_D_JAMES_128_Q5,
        fpoly: FP_5_94,
    },
    Params {
        name: "d-james-256-q13",
        scheme: Scheme::DJames,
        q: 13,
        n: 171,
        a: 53,
        m: 118,
        ny: 139,
        r: 2,
        d: 1,
        dd: 14,
        monomials: MON_D_JAMES_256_Q13,
        fpoly: FP_13_171,
    },
    Params {
        name: "d-james-256-q2",
        scheme: Scheme::DJames,
        q: 2,
        n: 390,
        a: 66,
        m: 324,
        ny: 512,
        r: 2,
        d: 2,
        dd: 5,
        monomials: MON_D_JAMES_256_Q2,
        fpoly: FP_2_390,
    },
    Params {
        name: "d-james-256-q23",
        scheme: Scheme::DJames,
        q: 23,
        n: 163,
        a: 53,
        m: 110,
        ny: 114,
        r: 2,
        d: 1,
        dd: 24,
        monomials: MON_D_JAMES_256_Q23,
        fpoly: FP_23_163,
    },
    Params {
        name: "d-james-256-q4",
        scheme: Scheme::DJames,
        q: 4,
        n: 223,
        a: 53,
        m: 170,
        ny: 256,
        r: 2,
        d: 2,
        dd: 17,
        monomials: MON_D_JAMES_256_Q4,
        fpoly: FP_4_223,
    },
    Params {
        name: "d-james-256-q5",
        scheme: Scheme::DJames,
        q: 5,
        n: 206,
        a: 53,
        m: 153,
        ny: 221,
        r: 2,
        d: 1,
        dd: 6,
        monomials: MON_D_JAMES_256_Q5,
        fpoly: FP_5_206,
    },
    Params {
        name: "james-128-q13",
        scheme: Scheme::James,
        q: 13,
        n: 91,
        a: 21,
        m: 70,
        ny: 0,
        r: 2,
        d: 1,
        dd: 14,
        monomials: MON_JAMES_128_Q13,
        fpoly: FP_13_91,
    },
    Params {
        name: "james-128-q2",
        scheme: Scheme::James,
        q: 2,
        n: 283,
        a: 27,
        m: 256,
        ny: 0,
        r: 2,
        d: 2,
        dd: 5,
        monomials: MON_JAMES_128_Q2,
        fpoly: FP_2_283,
    },
    Params {
        name: "james-128-q23",
        scheme: Scheme::James,
        q: 23,
        n: 78,
        a: 21,
        m: 57,
        ny: 0,
        r: 2,
        d: 1,
        dd: 24,
        monomials: MON_JAMES_128_Q23,
        fpoly: FP_23_78,
    },
    Params {
        name: "james-128-q4",
        scheme: Scheme::James,
        q: 4,
        n: 149,
        a: 21,
        m: 128,
        ny: 0,
        r: 2,
        d: 2,
        dd: 17,
        monomials: MON_JAMES_128_Q4,
        fpoly: FP_4_149,
    },
    Params {
        name: "james-128-q5",
        scheme: Scheme::James,
        q: 5,
        n: 132,
        a: 21,
        m: 111,
        ny: 0,
        r: 2,
        d: 1,
        dd: 6,
        monomials: MON_JAMES_128_Q5,
        fpoly: FP_5_132,
    },
    Params {
        name: "james-256-q13",
        scheme: Scheme::James,
        q: 13,
        n: 192,
        a: 53,
        m: 139,
        ny: 0,
        r: 2,
        d: 1,
        dd: 14,
        monomials: MON_JAMES_256_Q13,
        fpoly: FP_13_192,
    },
    Params {
        name: "james-256-q2",
        scheme: Scheme::James,
        q: 2,
        n: 578,
        a: 66,
        m: 512,
        ny: 0,
        r: 2,
        d: 2,
        dd: 5,
        monomials: MON_JAMES_256_Q2,
        fpoly: FP_2_578,
    },
    Params {
        name: "james-256-q23",
        scheme: Scheme::James,
        q: 23,
        n: 167,
        a: 53,
        m: 114,
        ny: 0,
        r: 2,
        d: 1,
        dd: 24,
        monomials: MON_JAMES_256_Q23,
        fpoly: FP_23_167,
    },
    Params {
        name: "james-256-q4",
        scheme: Scheme::James,
        q: 4,
        n: 309,
        a: 53,
        m: 256,
        ny: 0,
        r: 2,
        d: 2,
        dd: 17,
        monomials: MON_JAMES_256_Q4,
        fpoly: FP_4_309,
    },
    Params {
        name: "james-256-q5",
        scheme: Scheme::James,
        q: 5,
        n: 274,
        a: 53,
        m: 221,
        ny: 0,
        r: 2,
        d: 1,
        dd: 6,
        monomials: MON_JAMES_256_Q5,
        fpoly: FP_5_274,
    },
    Params {
        name: "toy-d-james-q13",
        scheme: Scheme::DJames,
        q: 13,
        n: 18,
        a: 4,
        m: 14,
        ny: 16,
        r: 1,
        d: 1,
        dd: 14,
        monomials: MON_TOY_D_JAMES_Q13,
        fpoly: FP_13_18,
    },
    Params {
        name: "toy-d-james-q2",
        scheme: Scheme::DJames,
        q: 2,
        n: 32,
        a: 6,
        m: 26,
        ny: 48,
        r: 2,
        d: 2,
        dd: 5,
        monomials: MON_TOY_D_JAMES_Q2,
        fpoly: FP_2_32,
    },
    Params {
        name: "toy-d-james-q23",
        scheme: Scheme::DJames,
        q: 23,
        n: 16,
        a: 4,
        m: 12,
        ny: 14,
        r: 1,
        d: 1,
        dd: 24,
        monomials: MON_TOY_D_JAMES_Q23,
        fpoly: FP_23_16,
    },
    Params {
        name: "toy-d-james-q4",
        scheme: Scheme::DJames,
        q: 4,
        n: 24,
        a: 5,
        m: 19,
        ny: 24,
        r: 2,
        d: 2,
        dd: 17,
        monomials: MON_TOY_D_JAMES_Q4,
        fpoly: FP_4_24,
    },
    Params {
        name: "toy-d-james-q5",
        scheme: Scheme::DJames,
        q: 5,
        n: 21,
        a: 4,
        m: 17,
        ny: 24,
        r: 2,
        d: 1,
        dd: 6,
        monomials: MON_TOY_D_JAMES_Q5,
        fpoly: FP_5_21,
    },
    Params {
        name: "toy-james-q13",
        scheme: Scheme::James,
        q: 13,
        n: 18,
        a: 4,
        m: 14,
        ny: 0,
        r: 1,
        d: 1,
        dd: 14,
        monomials: MON_TOY_JAMES_Q13,
        fpoly: FP_13_18,
    },
    Params {
        name: "toy-james-q2",
        scheme: Scheme::James,
        q: 2,
        n: 32,
        a: 6,
        m: 26,
        ny: 0,
        r: 2,
        d: 2,
        dd: 5,
        monomials: MON_TOY_JAMES_Q2,
        fpoly: FP_2_32,
    },
    Params {
        name: "toy-james-q23",
        scheme: Scheme::James,
        q: 23,
        n: 16,
        a: 4,
        m: 12,
        ny: 0,
        r: 1,
        d: 1,
        dd: 24,
        monomials: MON_TOY_JAMES_Q23,
        fpoly: FP_23_16,
    },
    Params {
        name: "toy-james-q4",
        scheme: Scheme::James,
        q: 4,
        n: 24,
        a: 5,
        m: 19,
        ny: 0,
        r: 2,
        d: 2,
        dd: 17,
        monomials: MON_TOY_JAMES_Q4,
        fpoly: FP_4_24,
    },
    Params {
        name: "toy-james-q5",
        scheme: Scheme::James,
        q: 5,
        n: 21,
        a: 4,
        m: 17,
        ny: 0,
        r: 2,
        d: 1,
        dd: 6,
        monomials: MON_TOY_JAMES_Q5,
        fpoly: FP_5_21,
    },
];

/// Look a parameter set up by name.
pub fn by_name(name: &str) -> Option<&'static Params> {
    ALL.iter().find(|p| p.name == name)
}

/// Names of every parameter set.
pub fn names() -> Vec<&'static str> {
    ALL.iter().map(|p| p.name).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn internally_consistent() {
        for p in ALL {
            assert_eq!(p.m, p.n - p.a, "{}: m = n - a", p.name);
            assert_eq!(p.q.pow(p.d as u32) + 1, p.dd as u32, "{}: D = q^d + 1", p.name);
            assert_eq!(p.fpoly.len(), p.n, "{}: field polynomial degree", p.name);
            assert!(p.fpoly[0] != 0, "{}: t divides f", p.name);
            assert!(p.fpoly.iter().all(|&c| (c as u32) < p.q), "{}: coefficient range", p.name);
            assert_eq!(p.is_dragon(), p.ny > 0, "{}: Dragon iff ny", p.name);
            for &(i, j) in p.monomials {
                assert!(i <= j && j <= p.d, "{}: monomial order", p.name);
                assert!(
                    p.q.pow(i as u32) + p.q.pow(j as u32) <= p.dd as u32,
                    "{}: monomial degree",
                    p.name
                );
            }
            // Over F_2, X^(q^i + q^i) is F_2-linear and useless as a quadratic
            // term, so no parameter set may use i = j there.
            if p.q == 2 {
                assert!(p.monomials.iter().all(|&(i, j)| i != j), "{}", p.name);
            }
        }
    }

    #[test]
    fn tags_are_well_formed() {
        assert_eq!(
            by_name("d-james-128-q2").unwrap().tag(),
            "D-James/v1/d-james/q2/n189/a27/r2/D5/ny256/mon0-1.0-2"
        );
        assert_eq!(
            by_name("james-128-q2").unwrap().tag(),
            "D-James/v1/james/q2/n283/a27/r2/D5/ny0/mon0-1.0-2"
        );
        assert_eq!(
            by_name("d-james-128-q23").unwrap().tag(),
            "D-James/v1/d-james/q23/n74/a21/r2/D24/ny57/mon0-0.0-1"
        );
        // Distinct sets must never share a tag: keys are bound to it.
        let mut tags: Vec<String> = ALL.iter().map(|p| p.tag()).collect();
        tags.sort();
        let n = tags.len();
        tags.dedup();
        assert_eq!(tags.len(), n, "duplicate parameter tag");
    }

    #[test]
    fn signature_sizes_match_the_spec() {
        // (name, bits, bytes) from d-james-spec.md section 9.2
        let want: &[(&str, usize)] = &[
            ("d-james-128-q2", 24), ("d-james-128-q4", 27), ("d-james-128-q5", 28),
            ("d-james-128-q13", 36), ("d-james-128-q23", 42), ("d-james-256-q2", 49),
            ("d-james-256-q4", 56), ("d-james-256-q5", 60), ("d-james-256-q13", 80),
            ("d-james-256-q23", 93), ("james-128-q2", 36), ("james-128-q4", 38),
            ("james-128-q5", 39), ("james-128-q13", 43), ("james-128-q23", 45),
            ("james-256-q2", 73), ("james-256-q4", 78), ("james-256-q5", 80),
            ("james-256-q13", 89), ("james-256-q23", 95),
        ];
        for &(name, bytes) in want {
            assert_eq!(by_name(name).unwrap().sig_bytes(), bytes, "{name}");
        }
    }
}

#!/usr/bin/env python3
"""Regenerate djames-rs/src/params.rs from the Python reference definitions.

    cd djames-py && python3 ../djames-rs/tools/gen_params_rs.py > ../djames-rs/src/params.rs
    cd ../djames-rs && cargo fmt

Generating rather than transcribing keeps the two implementations from drifting
apart through a typo in a 578-coefficient field polynomial. The file is emitted
whole, so nothing hand-written may live in it -- the parameter tests are in
tests/params.rs for exactly that reason.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "djames-py"))

from djames import params as PP  # noqa: E402

HEADER = '''//! Parameter sets and canonical field polynomials.
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

    /// The domain-separation tag of `d-james-spec.md` §2.4.
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
'''

FOOTER = '''
/// Look a parameter set up by name.
pub fn by_name(name: &str) -> Option<&'static Params> {
    ALL.iter().find(|p| p.name == name)
}

/// Names of every parameter set.
pub fn names() -> Vec<&'static str> {
    ALL.iter().map(|p| p.name).collect()
}'''


def main():
    fp = json.load(open(os.path.join(HERE, "..", "..", "djames-py",
                                     "djames", "fieldpoly.json")))
    names = PP.names()
    out = [HEADER]

    seen = {}
    for nm in names:
        p = PP.get(nm)
        seen.setdefault((p.q, p.n), fp["%d,%d" % (p.q, p.n)])
    for (q, n), cs in sorted(seen.items()):
        out.append("#[rustfmt::skip]\nstatic FP_%d_%d: &[u8] = &[%s];"
                   % (q, n, ", ".join(str(c) for c in cs)))
    out.append("")

    for nm in names:
        p = PP.get(nm)
        out.append("static MON_%s: &[(usize, usize)] = &[%s];"
                   % (nm.replace("-", "_").upper(),
                      ", ".join("(%d, %d)" % (i, j) for (i, j) in p.monomials)))
    out.append("")

    out.append("/// Every parameter set, in the order of `d-james-spec.md` §2.")
    out.append("pub static ALL: &[Params] = &[")
    order = ([n for n in names if not n.startswith("toy-")]
             + [n for n in names if n.startswith("toy-")])
    for nm in order:
        p = PP.get(nm)
        out.append("    Params { name: \"%s\", scheme: %s, q: %d, n: %d, a: %d,"
                   " m: %d, ny: %d, r: %d, d: %d, dd: %d, monomials: MON_%s,"
                   " fpoly: FP_%d_%d },"
                   % (nm,
                      "Scheme::DJames" if p.scheme == "d-james" else "Scheme::James",
                      p.q, p.n, p.a, p.m, p.ny or 0, p.r, p.d, p.D,
                      nm.replace("-", "_").upper(), p.q, p.n))
    out.append("];")
    out.append(FOOTER)
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()

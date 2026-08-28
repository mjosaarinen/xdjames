//! D-James and James: ultra-short multivariate signatures.
//!
//! An implementation of the schemes in Jacques Patarin and Alexandre Roullet,
//! *D-James: Ultra Short Multivariate Signatures*, IACR ePrint 2026/1650,
//! following `d-james-spec.md` in this repository. It is byte-for-byte
//! interoperable with the Python reference implementation in `djames-py/`:
//! both reproduce the same known-answer vectors.
//!
//! # Status
//!
//! This is a research implementation. The schemes themselves are new and the
//! paper's own authors recommend against deployment until they have had
//! substantially more public scrutiny.
//!
//! # Constant time
//!
//! Field arithmetic, the linear algebra over `F_q`, and public-key assembly
//! avoid secret-dependent branches and secret-indexed memory. Verification
//! handles only public data and is unconstrained. **Root finding is not
//! constant time** and cannot cheaply be made so: the number of roots of the
//! central polynomial and the control flow of equal-degree splitting depend on
//! secret-derived data. That is a property of HFE inversion rather than of
//! this implementation. See the crate README.
//!

#![cfg_attr(not(test), no_std)]
#![forbid(unsafe_op_in_unsafe_fn)]
// In indexed linear algebra the index *is* the meaning: `for i in 0..n` next
// to `for c in 0..n` mirrors the sum being implemented, and iterator adaptors
// obscure which axis is which. Kept deliberately.
#![allow(clippy::needless_range_loop)]

extern crate alloc;

pub mod codec;
pub mod ct;
pub mod gf;
pub mod keccak;
pub mod linalg;
pub mod params;
pub mod poly;
pub mod scheme;
pub mod symmetric;

pub use params::Params;
pub use scheme::{keygen, sign, verify, Error, PublicKey, SecretKey, SEED_BYTES};

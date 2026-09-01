# djames-rs

A performance-oriented, portable Rust implementation of **James** and
**D-James**, the ultra-short multivariate signature schemes of

> Jacques Patarin and Alexandre Roullet, *D-James: Ultra Short Multivariate
> Signatures*, IACR ePrint **2026/1650**, <https://eprint.iacr.org/2026/1650>

It follows [`../d-james-spec.md`](../d-james-spec.md) and is **byte-for-byte
interoperable** with the Python reference in [`../djames-py/`](../djames-py/):
`cargo test --release` re-derives that implementation's known-answer vectors —
public-key digests and every signature — for all five base fields and both
schemes.

> [!WARNING]
> Research code. The schemes are new and the paper's own authors "do not
> recommend deploying these signatures in security-critical applications until
> they have undergone further scrutiny." Root finding is not constant time
> (see below). Nothing here has had a security review.

## Using it

```rust
use djames::{params, keygen, sign, verify};

let p = params::by_name("d-james-128-q2").unwrap();
let (pk, sk) = keygen(p, &[1u8; 32])?;      // deterministic in the seed
let sig = sign(p, &sk, b"hello")?;          // 24 bytes = 189 bits
assert!(verify(p, &pk, b"hello", &sig));
```

```console
$ cargo test --release                  # unit tests + toy KAT vectors
$ cargo test --release -- --include-ignored  # + all 20 real parameter sets
$ cargo run --release --example bench
```

**No dependencies.** Not "few" — zero. SHAKE256, SHA3-256, the field
arithmetic and the big-integer base conversion are all in the crate, so it
builds anywhere with a 64-bit integer type. `#![no_std]` with `alloc`.

## Portability

There is no platform-specific code and no `unsafe` outside the two volatile
writes in `zeroize`. In particular the carry-less multiply that carries `F_2`
arithmetic is the portable four-way-nibble-split technique (BearSSL's
`bmul32`), not `pclmulqdq` or `vmull_p64`:

```rust
let (x0, x1, x2, x3) = (x & M0, x & M1, x & M2, x & M3);   // M0 = 0x1111...
let z0 = x0*y0 ^ x1*y3 ^ x2*y2 ^ x3*y1;                    // (four such)
```

Each output nibble accumulates at most eight terms, so an ordinary integer
multiply never carries between the classes kept. It is roughly 5–10× slower
than a hardware carry-less multiply; an intrinsic path would slot in behind a
`target_feature` gate at `gf::clmul64` and change nothing else.

## Representation

An element of `K = F_q^n` is a `Vec<u64>`, packed per field so the public key
stays in memory — `james-256-q2` holds 167 331 coefficient vectors, which is
11 MB packed against 86 MB at one byte per coefficient:

| q | packing | element |
|---|---|---|
| 2 | one coefficient per bit | `ceil(n/64)` words |
| 4 | two bit-planes | `2*ceil(n/64)` words |
| p odd | one coefficient per byte | `ceil(n/8)` words |

Reduction is nearly free. The canonical field polynomials are sparse with
every low exponent below 9, so writing `f = t^n + g`, the identity
`H = t^n G + L = gG + L` drops the degree from `2n-2` to at most `n+6` in one
fold and finishes in two. The fold count is fixed at construction, so
reduction is straight-line with no data-dependent indexing.

`F_4` products go through three `F_2[t]` products by Karatsuba, recombined
with `u² = u + 1`. Odd `F_p` products accumulate in `u32` without reduction —
the largest partial sum is `n(q-1)²`, far below `2^32` — then reduce once.

## Performance

The headline table is in the [root README](../README.md#performance);
`cargo run --release --example bench` prints fresh measurements and
`cargo run --release --example prof` breaks a signature down into its parts.

Two optimisations shape the signing cost:

**The Frobenius schedule is not one algorithm.** Computing `X^(q^n) mod F` by
doubling the Frobenius exponent costs `O(log n * D^3)`; iterating it costs
`O(n * D^2)`. Which wins flips with the parameter set — the ladder for
`q = 2` (`D = 5`, `n = 189`), the linear walk for `q = 13` and `q = 23`, where
`D` reaches 24 and `n` is only 74. `poly::x_pow_qn` selects per call.

**Equal-degree splitting over odd `F_q`** factors the exponent
`(q^n - 1)/2` as `(q-1)/2 * (1 + q + ... + q^(n-1))`, turning the inner part into
`prod_i b^(q^i)` — `n` Frobenius steps and `n` products, because raising to the
`q` inside `K[X]/(g)` is a linear combination against precomputed powers of
`X^q`, not a generic exponentiation.

The main `F_2` cost centre is the portable carry-less multiply, which is
where a hardware `pclmulqdq` / `vmull_p64` path would pay for itself.

## Constant time

The scheme divides cleanly into three parts.

**Public, so unconstrained.** Verification touches only the public key, the
signature and the message, so it skips zero coefficients freely.

**Constant time.** Base-field and extension-field arithmetic, the linear
algebra over `F_q` that handles `S`, `T` and `M_Z`, and public-key assembly.
No secret-dependent branch, no secret index, and no division by a secret —
`ct::barrett_u32` replaces `%` because integer division is variable-time on
several architectures. `F_q` has no lookup tables at all: a 529-entry table
indexed by secret values spans several cache lines and would be a timing
oracle, so `F_4` multiplication is bit arithmetic and prime fields use
Barrett. Elimination runs a fixed loop nest over every row and column, using
masks; for the full-rank square inputs `S` and `T` always are, the path is
independent of the entries. Secret key material is zeroized on drop.

**Not constant time, inherently: root finding.** The number of roots of the
central polynomial and the control flow of equal-degree splitting depend on
secret-derived data. That is a property of HFE inversion, not of this code —
no HFE implementation known to us avoids it. What is achievable, and done, is
that everything *underneath* it is constant time, so what leaks is the shape
of the factorisation rather than field-element values.

Two smaller items, named rather than hidden:

* Rejection sampling leaks only its number of rejections, which is independent
  of the accepted values; the matrix-rejection loops reject a *candidate*,
  never accepted key material.
* `random_full_rank` for `M_Z` leaves which columns became pivots dependent on
  the accepted matrix. For a uniform matrix these are the first `r` columns
  with overwhelming probability, so the leak is slight — but real.

`wrapping_mul` on small integers is assumed constant time, which holds on
mainstream 64-bit CPUs but is not a language guarantee; on a core with a
data-dependent multiplier, `gf::clmul64` would need a shift-and-mask
replacement. BearSSL's portable GHASH makes the same assumption.

## Layout

```
src/keccak.rs     Keccak-f[1600], SHAKE256, SHA3-256
src/symmetric.rs  domain-separated XOF, F_q sampling, message hashing
src/ct.rs         masks, selects, Barrett reduction, zeroize
src/gf.rs         F_q, K = F_q^n, portable carry-less multiply, packed vectors
src/poly.rs       polynomials over K; root finding (the trapdoor)
src/linalg.rs     matrices over F_q
src/params.rs     generated from djames-py; do not hand-edit
src/codec.rs      canonical encodings + the small bignum they need
src/scheme.rs     keygen / sign / verify
tests/kat.rs      the Python implementation's vectors
```

`src/params.rs` is generated by `tools/gen_params_rs.py` rather than
transcribed, so a 578-entry field polynomial cannot drift between the two
implementations through a typo.

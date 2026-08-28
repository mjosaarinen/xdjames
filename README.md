# xdjames

Two independent implementations and a specification for **James** and
**D-James**, the ultra-short multivariate signature schemes of

> Jacques Patarin and Alexandre Roullet, *D-James: Ultra Short Multivariate
> Signatures*, IACR Cryptology ePrint Archive **2026/1650**,
> <https://eprint.iacr.org/2026/1650>

D-James produces **24-byte signatures** at the 128-bit classical level — under
half of Ed25519's 64, and shorter than every post-quantum signature in the
paper's comparison — in exchange for a megabyte-scale public key. (The paper's
headline 156-bit figure is its truncated `|sig_short|` variant, whose
verification needs a hybrid Gröbner solve; that is not implemented here, so
24 bytes is the untruncated `|sig_fast|`.)

> [!WARNING]
> Research code, and a young scheme. The paper's own authors "do not recommend
> deploying these signatures in security-critical applications until they have
> undergone further scrutiny and extensive analysis by the cryptographic
> community." Nothing here has had a security review, and root finding is not
> constant time.

## Security target

The `*-128-*` parameter sets are intended to target **EUF-CMA security at
NIST security category 2**: an existential forgery under an adaptive
chosen-message attack should require resources comparable to collision search
on a 256-bit hash function.

This is a **provisional design claim**, not a proof or an independent security
assessment. It relies on the paper's estimates for hybrid, MinRank,
differential, and generic attacks; some of those estimates are heuristic and
the schemes have not yet received substantial third-party cryptanalysis. NIST's
baseline definition is EUF-CMA, with up to `2^64` chosen-message signing queries
in its estimation model; this repository currently makes no SUF-CMA claim.
Canonical decoding removes byte-encoding malleability, but odd-characteristic
James is algebraically not strongly unforgeable because negating a valid
nonzero signature produces a second valid signature for the same message.

See the [NIST PQC security criteria](https://csrc.nist.gov/Projects/Post-Quantum-Cryptography/Post-Quantum-Cryptography-Standardization/Evaluation-Criteria/Security-(Evaluation-Criteria))
and the warning above before interpreting the parameter names as established
security levels.

## What is here

| | |
|---|---|
| **[`d-james-spec.md`](d-james-spec.md)** | The specification. The paper gives the mathematics; this pins the hash, the field polynomials, the sampling order and the wire format — everything two implementations must agree on. |
| **[`djames-rs/`](djames-rs/)** | Rust. The one built for speed, portable and constant-time where achievable. Zero dependencies, `no_std`. |
| **[`djames-py/`](djames-py/)** | Python. The golden model: generates the test vectors and is meant to be read. Deliberately unoptimised, stdlib only. |
| **[`ref/`](ref/)** | Upstream material — the paper, the authors' SageMath proof-of-concept, and provenance for both. |

The two implementations were written against the specification and agree on
every shipped test vector, byte for byte. Coverage is uneven by design, because
signing cost explodes with `q` (see [Performance](#performance)):

| vectors | sets | what is pinned |
|---|---|---|
| `kat/toy.json` | 10 toy sets — all five fields, both schemes | public-key digest **and every signature** |
| `kat/q2.json` | the 4 real `F_2` sets | public-key digest **and every signature** |
| `kat/keygen.json` | all 20 real sets | public-key digest |

So full sign/verify agreement is established for 14 parameter sets and key
generation for all 30.

```console
$ make check          # both implementations, same vectors
$ make help           # everything else
```

Language-specific instructions live in
[`djames-rs/README.md`](djames-rs/README.md) and
[`djames-py/README.md`](djames-py/README.md).

## The schemes

**James** is HFE⁻_IP: a hidden low-degree univariate polynomial over `F_q^n`,
masked by two secret linear maps, with an *internal perturbation* to defeat the
MinRank-S attack that broke GeMSS and a *minus* modifier that publishes only
`m = n − a` of the `n` equations.

**D-James** adds **Dragon terms** — bilinear couplings `a_i b_j` between the
signature variables and the hash variables:

```
H(X, y) = Σ λ_ij X^(q^i + q^j)        HFE core
        + Σ L_k(y) · X^(q^k)          Dragon      ← D-James only
        + Σ M[k][j] · z_j · X^(q^k)   IP, bilinear
        + Σ G[i][j] · z_i z_j         IP, quadratic
```

Fixing the message hash `y` turns each `L_k(y)` into a *constant* of `K`, so
the central map collapses to a univariate polynomial of degree `D = q^d + 1`
and inverting it costs no more than plain HFE — while the number of published
equations `m` is decoupled from the hash length `n_y`. That decoupling is the
whole trick: it lets the signature be shorter than the `2λ` bits a hash-sized
signature would need.

## Parameter sets

`m = n − a` and `D = q^d + 1` throughout. Public-key sizes are the measured
serialized lengths.

### D-James — 24-byte signatures at 128-bit

| set | q | m | n | a | n_y | r | D | signature | public key |
|---|---|---|---|---|---|---|---|---|---|
| `d-james-128-q2` | 2 | 162 | 189 | 27 | 256 | 2 | 5 | **24 B** | 1,343,365 B |
| `d-james-128-q4` | 4 | 84 | 105 | 21 | 128 | 2 | 17 | 27 B | 399,105 B |
| `d-james-128-q5` | 5 | 73 | 94 | 21 | 111 | 2 | 6 | 28 B | 315,763 B |
| `d-james-128-q13` | 13 | 56 | 77 | 21 | 70 | 2 | 14 | 36 B | 218,218 B |
| `d-james-128-q23` | 23 | 53 | 74 | 21 | 57 | 2 | 24 | 42 B | 211,788 B |
| `d-james-256-q2` | 2 | 324 | 390 | 66 | 512 | 2 | 5 | 49 B | 11,174,963 B |
| `d-james-256-q4` | 4 | 170 | 223 | 53 | 256 | 2 | 17 | 56 B | 3,487,720 B |
| `d-james-256-q5` | 5 | 153 | 206 | 53 | 221 | 2 | 6 | 60 B | 2,969,301 B |
| `d-james-256-q13` | 13 | 118 | 171 | 53 | 139 | 2 | 14 | 80 B | 2,107,881 B |
| `d-james-256-q23` | 23 | 110 | 163 | 53 | 114 | 2 | 24 | 93 B | 2,008,160 B |

### James — no Dragon terms

| set | q | m | n | a | r | D | signature | public key |
|---|---|---|---|---|---|---|---|---|
| `james-128-q2` | 2 | 256 | 283 | 27 | 2 | 5 | 36 B | 1,285,952 B |
| `james-128-q4` | 4 | 128 | 149 | 21 | 2 | 17 | 38 B | 357,600 B |
| `james-128-q5` | 5 | 111 | 132 | 21 | 2 | 6 | 39 B | 282,879 B |
| `james-128-q13` | 13 | 70 | 91 | 21 | 2 | 14 | 43 B | 136,045 B |
| `james-128-q23` | 23 | 57 | 78 | 21 | 2 | 24 | 45 B | 100,353 B |
| `james-256-q2` | 2 | 512 | 578 | 66 | 2 | 5 | 73 B | 10,709,184 B |
| `james-256-q4` | 4 | 256 | 309 | 53 | 2 | 17 | 78 B | 3,065,280 B |
| `james-256-q5` | 5 | 221 | 274 | 53 | 2 | 6 | 80 B | 2,417,277 B |
| `james-256-q13` | 13 | 139 | 192 | 53 | 2 | 14 | 89 B | 1,195,718 B |
| `james-256-q23` | 23 | 114 | 167 | 53 | 2 | 24 | 95 B | 913,824 B |

Ten further *toy* sets (`toy-*`) exist for tests only — no security whatsoever;
see [spec §2.3](d-james-spec.md).

Every signature size matches the paper's `|sig_fast|` column exactly, bar the
`q = 23` rows where the paper reports 334 and 737 bits against `⌈n log₂ 23⌉` of
335 and 738. The public-key column is a different story: it reproduces the
paper's James tables to within a few percent but not its D-James ones — see
[`djames-py/README.md`](djames-py/README.md#discrepancies-in-the-papers-tables).

## Performance

[`djames-rs/`](djames-rs/), single-threaded, on a 24-core x86-64 box. Sizes are
measured serialized lengths. `cargo run --release --example bench`.

| set | q | n | keygen | sign | verify | signature | public key |
|---|---|---|---|---|---|---|---|
| `d-james-128-q2` | 2 | 189 | 244 ms | 93 ms | **131 µs** | **24 B** | 1.34 MB |
| `james-128-q2` | 2 | 283 | 332 ms | 226 ms | 58 µs | 36 B | 1.29 MB |
| `d-james-256-q2` | 2 | 390 | 1.4 s | 1.2 s | 447 µs | 49 B | 11.2 MB |
| `james-256-q2` | 2 | 578 | 2.6 s | 2.8 s | 321 µs | 73 B | 10.7 MB |
| `d-james-128-q4` | 4 | 105 | 603 ms | 1.1 s | 3.2 ms | 27 B | 399 kB |
| `d-james-128-q5` | 5 | 94 | 704 ms | 1.2 s | 4.5 ms | 28 B | 316 kB |
| `d-james-128-q13` | 13 | 77 | 302 ms | 2.9 s | 2.7 ms | 36 B | 218 kB |
| `d-james-128-q23` | 23 | 74 | 254 ms | 244 s | 7.4 ms | 42 B | 212 kB |

Two things that column is telling you.

**Verification is cheap and signing is not, and the gap widens with `q`.** The
IP modifier costs `q^r` root-findings per salt attempt, and every parameter set
uses `r = 2` — so `q = 2` pays a factor of 4 while `q = 23` pays **529**, each
one a degree-24 root-finding over `F_{23^74}`. That is the scheme, not the
implementation: the paper says outright that "the parameter `q^r` must remain
small in practice", and its own performance estimates (Table 9) cover the `F_2`
sets alone. If you want D-James to be usable, `q = 2` is the row to look at.

**The `F_2` sets are the practical ones anyway** — they carry the headline
24-byte signature, and they verify in about 130 µs against a megabyte-scale
key.

Two caveats worth stating. The carry-less multiply is the portable
nibble-split construction rather than `pclmulqdq`, which costs roughly 5–10× on
the `F_2` sets; an intrinsic path would slot in behind a `target_feature` gate
without touching anything else. And these are single-sample timings on a busy
machine, not a careful benchmark — treat them as orders of magnitude.


## File map

```
d-james-spec.md          the specification: parameters, fields, pseudocode, wire format
Makefile                 make check / test / kat / bench / clean  (make help)

djames-rs/               Rust, performance and constant time
  src/keccak.rs            Keccak-f[1600], SHAKE256, SHA3-256
  src/symmetric.rs         domain-separated XOF, F_q sampling, message hashing
  src/ct.rs                masks, selects, Barrett reduction, zeroize
  src/gf.rs                F_q, K = F_q^n, portable carry-less multiply, packed vectors
  src/poly.rs              polynomials over K; root finding (the trapdoor)
  src/linalg.rs            matrices over F_q
  src/codec.rs             canonical encodings + the small bignum they need
  src/params.rs            generated from djames-py; do not hand-edit
  src/scheme.rs            keygen / sign / verify
  tests/kat.rs             the golden model's vectors
  examples/bench.rs        timings

djames-py/               Python, the golden model
  djames/ff.py             F_q, packed backends, K = F_q^n, irreducibility
  djames/poly.py           polynomials over K; root finding
  djames/linalg.py         matrices over F_q
  djames/symmetric.py      SHAKE256 XOF, F_q sampling
  djames/codec.py          canonical encodings
  djames/scheme.py         keygen / sign / verify
  djames/params.py         parameter sets
  djames/fieldpoly.json    cached canonical field polynomials
  djames/kat.py            vector generation and checking
  kat/*.json, *.rsp        the test vectors, in both formats
  tools/                   regenerate vectors, field polynomials, params.rs

ref/                     upstream: the paper, the authors' Sage notebook, provenance
```

## Reading the paper

Three of the paper's statements are self-contradictory, and the specification
resolves each in favour of the reading the rest of the paper supports. All
three would benefit from the authors' confirmation:

1. **Dragon index range** — the text names `L_0, …, L_{d−1}` but the displayed
   sum runs over `q^i ≤ D`, which includes `i = d`.
2. **`q = 4` central monomials** — all four `q = 4` rows give `X³, X¹⁷`, but
   over `F_4` there is no `X³ = X^(q^i + q^j)`.
3. **`m` for `q = 4` at 256-bit** — Table 4 lists `m = 180` with `n = 223`,
   `a = 53`, and Table 8 repeats it, against the header's own `m = n − a`.

Details in [spec §11](d-james-spec.md).

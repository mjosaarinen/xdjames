# `ref/` — upstream reference material

Provenance for everything in this directory. Nothing here is our own work
except `djames_demo.sage`, which is a mechanical extraction (see below).

## `HFE_companion_notebook.ipynb`

The authors' proof-of-concept, cited in §7.1 of the paper.

- Source: <https://github.com/alex-on-the-internet/D-James-Sage-implementation>
- Retrieved from branch `main`, 2026-08-28
- Committed by `alexandreroullet3-dotcom` (the paper's second author),
  latest commit 2026-08-12T13:06:38Z ("Add files via upload"); 3 commits total
- Repository description: *"This is a proof-of-concept implementation of
  D-James and James, it has no intention for secure cryptographic
  applications."*
- **No license file** in the upstream repository — treat as all-rights-reserved
  and do not copy code from it into anything we ship.
- `sha256: 82c6722ebbbcf6e5572519e6737a33af85e7e646b05bd468d3c8ef7aa00358b5`
  (byte-identical to upstream `main`)

SageMath 10.7 kernel, 25 cells. Implements four schemes — plain HFE,
HFE-Dragon, HFE-IP, and HFE-IP-Dragon (= D-James) — as `keygen`/`sign`/`verify`
triples, plus a helper that exports the public key as explicit multivariate
quadratic polynomials over `GF(q)`.

What it is *not*: the demo parameters are toys (`n=48, m=32, ny=64` for q=2 and
`n=21, m=14, ny=28` for q=5) and none of the paper's parameter sets are
exercised; the hash is a bare SHA-256 digit expansion with no domain
separation; nothing is constant-time. The differential-attack code behind the
paper's Table 1 and Figure 1 is **not** published anywhere.

## `djames_demo.sage`

Ours, but derived: the notebook's code cells flattened into a single runnable
Sage script. The only edits are mechanical — Jupyter magics (`%display latex`,
`reset()`) stripped, and the two demo cells (22, 23) replaced by direct calls
that exercise HFE-Dragon and D-James at the notebook's own toy parameters.
No algorithmic change.

- `sha256: 840395e723db0e4d997f62c4296dad167429745b2822ca71e98fc9297b2eefcf`
- Run with: `sage djames_demo.sage`

Confirmed working on 2026-08-28 (SageMath 10.7):

```
--- HFE-Dragon (q=2)  (q=2, n=48, m=32, hash_length=64) ---
  trial 1: verify(sign(message)) = True
  trial 2: verify(sign(message)) = True
  tampered signature accepted = False  (expected: False)

--- D-James=HFE-IP-Dragon (q=2)  (q=2, n=48, m=32, r=2, hash_length=64) ---
  trial 1: verify(sign(message)) = True
  trial 2: verify(sign(message)) = True
  tampered signature accepted = False  (expected: False)

--- D-James=HFE-IP-Dragon (q=5)  (q=5, n=21, m=14, r=2, hash_length=28) ---
  trial 1: verify(sign(message)) = True
  trial 2: verify(sign(message)) = True
  tampered signature accepted = False  (expected: False)
```

## Relationship to `../djames-py/`

`djames-py` is an independent clean-room implementation written from the paper
and the normative repository specification, not from the notebook. It is
pure-Python, stdlib-only, uses SHAKE256 for hashing and pseudorandomness, and
produces deterministic test vectors.

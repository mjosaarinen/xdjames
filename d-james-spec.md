# D-James and James — Algorithm Specification

**Version 1.** This document specifies the two signature schemes of

> Jacques Patarin and Alexandre Roullet, *D-James: Ultra Short Multivariate
> Signatures*, IACR ePrint **2026/1650**, <https://eprint.iacr.org/2026/1650>

in enough detail to write an interoperable implementation. The paper defines
the mathematics; it does not fix a hash, a field polynomial, a sampling order,
or a wire format. Everything an implementation must agree on is pinned here.

Two independent implementations exist and agree on every test vector:
[`djames-py/`](djames-py/) (reference, readability first) and
[`djames-rs/`](djames-rs/) (performance, constant-time where achievable).

> [!WARNING]
> Neither the schemes nor this specification have had meaningful public
> review. The paper's authors "do not recommend deploying these signatures in
> security-critical applications until they have undergone further scrutiny."
> Treat this as a specification for study, not for deployment.

---

## 1. Notation

| symbol | meaning |
|---|---|
| `q` | base field size, one of 2, 4, 5, 13, 23 |
| `F_q` | the base field; elements are integers in `[0, q)` |
| `n` | degree of the extension `K = F_q^n`; also the signature length in F_q symbols |
| `a` | number of *minus* modifiers |
| `m` | number of published equations, `m = n − a` |
| `r` | rank of the *internal perturbation* (IP) modifier |
| `D` | degree of the central polynomial, `D = q^d + 1` |
| `d` | the exponent above |
| `n_y` | hash length in F_q symbols (D-James only) |
| `α` | a fixed generator of `K` over `F_q`; `{1, α, …, α^(n−1)}` is the working basis |
| `‖` | byte concatenation |
| `u32le(x)`, `u64le(x)` | `x` as 4 resp. 8 little-endian bytes |

Vectors are **rows**. For a matrix `M`, `v·M` denotes the row-vector product
`(v·M)_j = Σ_i v_i M[i][j]`, and `M·v` the column product
`(M·v)_i = Σ_j M[i][j] v_j`.

An element of `K` is identified with its coordinate vector over the basis
above: `X = Σ_{i<n} x_i α^i` ↔ `x = (x_0, …, x_{n−1}) ∈ F_q^n`.

---

## 2. Parameters

`m = n − a` always. `d` is determined by `D = q^d + 1`. The `monomials` column
lists the pairs `(i, j)` with a nonzero `λ_{i,j} X^(q^i + q^j)` in the central
map, given both as index pairs and as the resulting exponent.

### 2.1 D-James

| set | q | m | n | a | n_y | r | D | d | monomials |
|---|---|---|---|---|---|---|---|---|---|
| `d-james-128-q2` | 2 | 162 | 189 | 27 | 256 | 2 | 5 | 2 | (0,1),(0,2) = X³, X⁵ |
| `d-james-128-q4` | 4 | 84 | 105 | 21 | 128 | 2 | 17 | 2 | (0,0),(0,2) = X², X¹⁷ |
| `d-james-128-q5` | 5 | 73 | 94 | 21 | 111 | 2 | 6 | 1 | (0,0),(0,1) = X², X⁶ |
| `d-james-128-q13` | 13 | 56 | 77 | 21 | 70 | 2 | 14 | 1 | (0,0),(0,1) = X², X¹⁴ |
| `d-james-128-q23` | 23 | 53 | 74 | 21 | 57 | 2 | 24 | 1 | (0,0),(0,1) = X², X²⁴ |
| `d-james-256-q2` | 2 | 324 | 390 | 66 | 512 | 2 | 5 | 2 | (0,1),(0,2) |
| `d-james-256-q4` | 4 | 170 | 223 | 53 | 256 | 2 | 17 | 2 | (0,0),(0,2) |
| `d-james-256-q5` | 5 | 153 | 206 | 53 | 221 | 2 | 6 | 1 | (0,0),(0,1) |
| `d-james-256-q13` | 13 | 118 | 171 | 53 | 139 | 2 | 14 | 1 | (0,0),(0,1) |
| `d-james-256-q23` | 23 | 110 | 163 | 53 | 114 | 2 | 24 | 1 | (0,0),(0,1) |

### 2.2 James

| set | q | m | n | a | r | D | d | monomials |
|---|---|---|---|---|---|---|---|---|
| `james-128-q2` | 2 | 256 | 283 | 27 | 2 | 5 | 2 | (0,1),(0,2) |
| `james-128-q4` | 4 | 128 | 149 | 21 | 2 | 17 | 2 | (0,0),(0,2) |
| `james-128-q5` | 5 | 111 | 132 | 21 | 2 | 6 | 1 | (0,0),(0,1) |
| `james-128-q13` | 13 | 70 | 91 | 21 | 2 | 14 | 1 | (0,0),(0,1) |
| `james-128-q23` | 23 | 57 | 78 | 21 | 2 | 24 | 1 | (0,0),(0,1) |
| `james-256-q2` | 2 | 512 | 578 | 66 | 2 | 5 | 2 | (0,1),(0,2) |
| `james-256-q4` | 4 | 256 | 309 | 53 | 2 | 17 | 2 | (0,0),(0,2) |
| `james-256-q5` | 5 | 221 | 274 | 53 | 2 | 6 | 1 | (0,0),(0,1) |
| `james-256-q13` | 13 | 139 | 192 | 53 | 2 | 14 | 1 | (0,0),(0,1) |
| `james-256-q23` | 23 | 114 | 167 | 53 | 2 | 24 | 1 | (0,0),(0,1) |

See §11 for the three places where these differ from the paper's tables.

### 2.3 Toy sets

**No security whatsoever.** These exist only so a test run can exercise every
field and both schemes in seconds; they are listed because the shipped test
vectors reference them by name. Note `r = 1` for `q = 13` and `q = 23`: signing
costs `q^r` root-findings and 529 of them is not a unit test.

| set | q | m | n | a | n_y | r | D | d | monomials |
|---|---|---|---|---|---|---|---|---|---|
| `toy-d-james-q2` | 2 | 26 | 32 | 6 | 48 | 2 | 5 | 2 | (0,1),(0,2) |
| `toy-d-james-q4` | 4 | 19 | 24 | 5 | 24 | 2 | 17 | 2 | (0,0),(0,2) |
| `toy-d-james-q5` | 5 | 17 | 21 | 4 | 24 | 2 | 6 | 1 | (0,0),(0,1) |
| `toy-d-james-q13` | 13 | 14 | 18 | 4 | 16 | 1 | 14 | 1 | (0,0),(0,1) |
| `toy-d-james-q23` | 23 | 12 | 16 | 4 | 14 | 1 | 24 | 1 | (0,0),(0,1) |
| `toy-james-q2` | 2 | 26 | 32 | 6 | — | 2 | 5 | 2 | (0,1),(0,2) |
| `toy-james-q4` | 4 | 19 | 24 | 5 | — | 2 | 17 | 2 | (0,0),(0,2) |
| `toy-james-q5` | 5 | 17 | 21 | 4 | — | 2 | 6 | 1 | (0,0),(0,1) |
| `toy-james-q13` | 13 | 14 | 18 | 4 | — | 1 | 14 | 1 | (0,0),(0,1) |
| `toy-james-q23` | 23 | 12 | 16 | 4 | — | 1 | 24 | 1 | (0,0),(0,1) |

Their field polynomials, by the rule of §3.2:

| q | n | nonzero c_i |
|---|---|---|
| 2 | 32 | c₀=1, c₂=1, c₃=1, c₇=1 |
| 4 | 24 | c₀=2, c₁=1, c₃=2, c₄=1 |
| 5 | 21 | c₀=1, c₁=4 |
| 13 | 18 | c₀=2 |
| 23 | 16 | c₀=20, c₁=1 |

### 2.4 Parameter tag

Key derivation is bound to the full parameter set through a printable-ASCII
tag with an explicit grammar (no language's `repr`):

```
D-James/v1/<scheme>/q<q>/n<n>/a<a>/r<r>/D<D>/ny<ny>/mon<i>-<j>[.<i>-<j>]...
```

`<scheme>` is `d-james` or `james`; every number is plain decimal with no
padding and no spaces; `ny` is `0` for James; monomial pairs are separated by
`.` in the order listed above. Examples:

```
D-James/v1/d-james/q2/n189/a27/r2/D5/ny256/mon0-1.0-2
D-James/v1/james/q2/n283/a27/r2/D5/ny0/mon0-1.0-2
D-James/v1/d-james/q23/n74/a21/r2/D24/ny57/mon0-0.0-1
```

---

## 3. Fields

### 3.1 The base field `F_q`

For prime `q` (2, 5, 13, 23), `F_q = Z/qZ` and an element is its residue in
`[0, q)`.

For `q = 4`, `F_4 = F_2[u]/(u² + u + 1)`. An element `e ∈ [0, 4)` encodes
`e = b_0 + 2b_1` ↔ `b_0 + b_1 u`. Addition is XOR of the two bits;
multiplication follows from `u² = u + 1`.

Generally, for `q = p^k` an element encodes its base-`p` digit vector as the
coefficients of a polynomial in `u`, reduced modulo the lexicographically
least monic irreducible of degree `k` over `F_p`. For `k = 2, p = 2` that is
`u² + u + 1`.

### 3.2 The extension field `K = F_q^n`

`K = F_q[t] / (f(t))` with `f` monic irreducible of degree `n`, and `α` is the
class of `t`.

**Choice of `f`.** Write `f = t^n + Σ_{i<n} c_i t^i`. Enumerate candidates by
the integer whose base-`q` digits are `(c_0, …, c_{n−1})`, starting at 1,
skipping any with `c_0 = 0` (those are divisible by `t`); `f` is the first
irreducible candidate. The rule is deterministic, so the table below is a
convenience, not an extra input.

| q | n | nonzero coefficients of `f = t^n + Σ c_i t^i` |
|---|---|---|
| 2 | 189 | c₀=1, c₂=1, c₅=1, c₆=1 |
| 2 | 283 | c₀=1, c₁=1, c₂=1, c₅=1, c₆=1, c₈=1 |
| 2 | 390 | c₀=1, c₁=1, c₃=1, c₄=1, c₅=1, c₈=1 |
| 2 | 578 | c₀=1, c₁=1, c₂=1, c₅=1, c₆=1, c₈=1 |
| 4 | 105 | c₀=1, c₄=1 |
| 4 | 149 | c₀=2, c₁=1, c₃=1 |
| 4 | 223 | c₀=3, c₂=3, c₃=1, c₄=2 |
| 4 | 309 | c₀=2, c₂=2, c₃=1, c₄=1 |
| 5 | 94 | c₀=1, c₁=1 |
| 5 | 132 | c₀=1, c₂=3, c₃=1 |
| 5 | 206 | c₀=3, c₂=1, c₃=1 |
| 5 | 274 | c₀=1, c₁=1, c₂=4, c₃=1 |
| 13 | 77 | c₀=2, c₁=4 |
| 13 | 91 | c₀=8, c₂=1 |
| 13 | 171 | c₀=3, c₁=1, c₂=1 |
| 13 | 192 | c₀=2 |
| 23 | 74 | c₀=3, c₁=9, c₂=1 |
| 23 | 78 | c₀=3, c₁=3, c₂=1 |
| 23 | 163 | c₀=3, c₁=1 |
| 23 | 167 | c₀=7, c₁=5, c₂=4 |

**Frobenius.** Since every coordinate lies in `F_q` and is fixed by `x ↦ x^q`,

```
(Σ_i a_i α^i)^(q^k) = Σ_i a_i (α^i)^(q^k) = Σ_i a_i (α^(q^k))^i
```

so raising to the `q^k` is an `F_q`-linear combination against the precomputed
table `[(α^(q^k))^i]_{i<n}`.

---

## 4. Symmetric layer

SHAKE256 is the only symmetric primitive. Every stream is opened under one of
five domain labels:

```
DOM_KEY = "D-James/v1/keygen"     expand the master seed
DOM_PRF = "D-James/v1/prf"        derive the per-key signing PRF key
DOM_MSG = "D-James/v1/msg"        hash a message to F_q^{n_y} (or F_q^m)
DOM_PAD = "D-James/v1/pad"        James: the random minus-padding
DOM_EDF = "D-James/v1/edf"        randomness inside equal-degree splitting
```

### 4.1 XOF

```
XOF(part_0, part_1, …):
    return SHAKE256 stream over  ‖_i ( u32le(len(part_i)) ‖ part_i )
```

Each part is length-prefixed so that `("ab","c")` and `("a","bc")` are distinct
inputs. `XOF.read(k)` returns the next `k` bytes of the squeezed stream.

### 4.2 Sampling `F_q`

```
SampleFq(xof, count, q) -> [F_q; count]:
    if q is a power of two:
        k    <- log2(q)
        buf  <- xof.read( ceil(count*k / 8) )          # exact
        acc  <- little-endian integer of buf
        for i in 0..count-1:
            out[i] <- acc & (q-1);  acc <- acc >> k
    else:
        limit <- floor(256/q) * q
        while len(out) < count:
            b <- xof.read(1)[0]                        # ONE byte at a time
            if b < limit:  out.append(b mod q)
    return out
```

Both paths consume a precisely determined number of stream bytes, which
matters because callers keep drawing from the same XOF afterwards. Reading the
rejection path in larger chunks yields the same digits but advances the stream
by a chunk-size-dependent amount, so **byte-at-a-time is normative.**

```
SampleK(xof)    -> K:  from_coords( SampleFq(xof, n, q) )
SampleKNonzero  -> K:  repeat SampleK until the result is nonzero
```

### 4.3 Hashing a message

```
HashToFq(msg, salt, count, q):
    x <- XOF(DOM_MSG, msg, u64le(salt), u32le(q), u32le(count))
    return SampleFq(x, count, q)
```

---

## 5. The central map

With `z = M_Z · x ∈ F_q^r` and `L_k(y) = Σ_{j<n_y} Λ[k][j] · y_j ∈ K`:

```
H(X, y) = Σ_{(i,j) ∈ monomials}  λ_{i,j} X^(q^i + q^j)          HFE core
        + Σ_{k = 0..d}           L_k(y) · X^(q^k)               Dragon  [D-James only]
        + Σ_{k = 0..d-1} Σ_{j<r} M[k][j] · z_j · X^(q^k)        IP, bilinear
        + Σ_{0 ≤ i ≤ j < r}      G[i][j] · z_i · z_j            IP, quadratic
```

Every term is quadratic in `(a, b)` jointly, so **the public system is
homogeneous**: it has no constant and no linear part, and the all-zero vector
is always a solution (see §7).

Fixing `y` turns each `L_k(y)` into a constant of `K`, so `H` collapses to a
univariate polynomial of degree `D` and inverting it costs no more than plain
HFE — while `m` is decoupled from `n_y`. That is the point of the construction:
the signature can be shorter than `2λ` bits.

---

## 6. Key generation

The secret key **is** its 32-byte seed; everything below is derived. Seeds of
any other length MUST be rejected.

```
KeyGen(P, seed[32]) -> (pk, sk):
    xof <- XOF(DOM_KEY, seed, P.tag())

    # --- draws, in exactly this order ---
    for t in 0 .. len(monomials)-1:   λ[t]       <- SampleKNonzero(xof)
    for k in 0 .. d-1:  for j in 0..r-1:   M[k][j]  <- SampleK(xof)
    for i in 0 .. r-1:  for j in i..r-1:   G[i][j]  <- SampleK(xof)
    if D-James:
        for k in 0 .. d:  for j in 0..n_y-1:  Λ[k][j] <- SampleK(xof)
    M_Z        <- RandomFullRank(xof, r, n)
    (S, S⁻¹)   <- RandomInvertible(xof, n)
    (T, T⁻¹)   <- RandomInvertible(xof, n)

    prf <- XOF(DOM_PRF, seed, P.tag()).read(32)
    pk  <- BuildPublicKey(...)
```

`G[i][j]` is drawn only for `i ≤ j`; the loop visits `(0,0),(0,1),…,(0,r−1),
(1,1),…` and consumes nothing for `i > j`.

```
RandomMatrix(xof, rows, cols):     # row-major
    flat <- SampleFq(xof, rows*cols, q);  reshape to rows x cols

RandomInvertible(xof, n):
    loop:  M <- RandomMatrix(xof, n, n)
           if M is invertible: return (M, M⁻¹)

RandomFullRank(xof, rows, cols):
    loop:  M <- RandomMatrix(xof, rows, cols)
           if rank(M) = rows: return M
```

Rejection is over the *sampled* matrix, so the number of retries is
independent of the accepted value. The inverse and the rank are unique, so the
elimination algorithm is unconstrained.

### 6.1 Building the public key

Write the central map's inputs as linear forms in `(a, b)` with `K`
coefficients. Since `x = a·S`:

```
σ_i      = Σ_c S[i][c] · α^c                       # X = Σ_i a_i σ_i
A[k][i]  = σ_i^(q^k)            for k = 0..d       # X^(q^k) = Σ_i a_i A[k][i]
ζ[j][i]  = Σ_c M_Z[j][c] · S[i][c]  ∈ F_q          # z_j    = Σ_i a_i ζ[j][i]
```

Each piece of `H` is then a product of two linear forms, i.e. a rank-one outer
product `u ⊗ v` over `K`. Folding the IP families over their inner index first
cuts the term count from `d·r + r(r+1)/2` to `d + r`:

```
aa_terms:                                                  # both sides in a
    for each monomial (i,j) with coefficient λ:
        ( [λ · A[i][p]]_p ,  A[j] )
    for k in 0..d-1:
        ( A[k] ,  [Σ_j M[k][j] · ζ[j][p]]_p )
    for i in 0..r-1:
        ( [ζ[i][p]]_p ,  [Σ_{j≥i} G[i][j] · ζ[j][p]]_p )

ab_terms:                                                  # a against b
    for k in 0..d:   ( A[k] , Λ[k] )                        [D-James only]
```

Accumulate the `K`-valued coefficient of each monomial, then project once:

```
C_aa[p][s] = Σ_terms ( u[p]·v[s] + u[s]·v[p] )   for p < s
C_aa[p][p] = Σ_terms   u[p]·v[p]
C_ab[p][s] = Σ_terms   u[p]·v[s]

π(w) = ( T·coords(w) )[0 .. m-1]        # projection + minus modifier
pk.aa[p][s] = π(C_aa[p][s]);   pk.ab[p][s] = π(C_ab[p][s])
```

The public key is `m` homogeneous quadratic equations over `F_q`:

```
Σ_{p ≤ s} pk.aa[p][s]·a_p a_s  +  Σ_{p,s} pk.ab[p][s]·a_p b_s  =  0
```

---

## 7. Signing

Deterministic: no entropy beyond the key and the message.

```
Sign(P, sk, msg) -> sig:
    edf <- XOF(DOM_EDF, sk.prf, msg)

    for counter = 0, 1, 2, … , MAX_SALT-1:

        if D-James:
            y      <- HashToFq(msg, counter, n_y, q)
            c_k    <- Σ_{j<n_y} Λ[k][j]·y_j        for k = 0..d
            target <- 0                             # solve H = 0
        else:                                       # James
            h      <- HashToFq(msg, 0, m, q)        # salt fixed at 0
            pad    <- SampleFq( XOF(DOM_PAD, sk.prf, msg, u64le(counter)), a, q )
            u      <- T⁻¹ · (h ‖ pad)
            target <- Σ_i u_i α^i                   # solve H = target

        for e = 0 .. q^r - 1:                       # the IP guess
            z_j <- floor(e / q^(r-1-j)) mod q       for j = 0..r-1
            F   <- CentralPoly(λ, Λ or c_k, M, G, z) - target
            for X in Roots(F, edf):                 # canonical order, §7.1
                if X = 0: continue                  # homogeneous: trivial root
                x <- coords(X)
                if M_Z·x ≠ z: continue              # IP guess was wrong
                A <- x·S⁻¹
                if A ≠ 0: return Encode(A)
    fail
```

`MAX_SALT = 256`. Each salt succeeds with probability about `1 − 1/e`
(characteristic 2) or `1 − 1/3.07` (odd), so the expected count is the paper's
1.582 resp. 3.07 and 256 salts fail with probability about `2⁻³⁶⁹` resp.
`2⁻¹⁴⁶`. The bound is **security-relevant**: a forger may aim at any of
`MAX_SALT` hash values, conceding `log2(MAX_SALT) = 8` bits.

`CentralPoly` places coefficients at the exponents named in §5: `λ_{i,j}` at
`q^i + q^j`, `c_k` and the IP bilinear sums at `q^k`, and the IP quadratic sum
`Σ_{i≤j} G[i][j] z_i z_j` at exponent 0.

### 7.1 Root ordering is normative

`Roots(F)` returns **all** roots of `F` in `K`, sorted by their coordinate
vectors read as base-`q` integers with `c_{n−1}` most significant.

Any correct root-finding algorithm may be used — typically
`gcd(X^(q^n) − X, F)` followed by equal-degree splitting. But equal-degree
splitting returns roots in an order that depends on the random `δ` it happens
to draw, and the signer takes the *first* root passing the IP check. Without a
canonical order two conforming implementations could emit different — both
valid — signatures for the same key and message. The `DOM_EDF` stream is used
only inside root finding, so how much of it an implementation consumes does
not affect anything else.

---

## 8. Verification

```
Verify(P, pk, msg, sig) -> bool:
    A <- Decode(sig)               # rejects wrong length or non-canonical, §9
    if A = 0: return false         # the trivial root of a homogeneous system

    if D-James:
        for salt = 0 .. MAX_SALT-1:
            b <- HashToFq(msg, salt, n_y, q)
            if Eval(pk, A, b) = 0: return true
        return false
    else:
        h <- HashToFq(msg, 0, m, q)
        return Eval(pk, A, ⊥) = h
```

The salt is not transmitted (the paper's footnote 5 trades verification time
for a shorter signature), so the verifier re-derives it by trying
`0, 1, 2, …`. The signer scans salts in the same order, so the expected cost
is under two evaluations.

---

## 9. Serialization

All wire formats are little-endian and canonical: **a decoder MUST reject any
byte string that is not the unique encoding of its value.**

### 9.1 Digit vectors

Let `vec_bytes(count, q) = ceil(count · log2 q / 8)`, computed exactly as
`ceil(bitlen(q^count − 1) / 8)`.

```
EncodeVector(digits, q):
    v <- Σ_i digits[i] · q^i
    return v as vec_bytes(len(digits), q) little-endian bytes

DecodeVector(data, count, q):
    require len(data) = vec_bytes(count, q)
    v <- little-endian integer of data
    require v < q^count                      # rejects non-canonical encodings
    return [ (v / q^i) mod q ]_{i<count}
```

### 9.2 Signature

`EncodeVector(A, q)` over the `n` symbols — exactly `vec_bytes(n, q)` bytes.

| set | bits | bytes | | set | bits | bytes |
|---|---|---|---|---|---|---|
| `d-james-128-q2` | 189 | **24** | | `james-128-q2` | 283 | 36 |
| `d-james-128-q4` | 210 | 27 | | `james-128-q4` | 298 | 38 |
| `d-james-128-q5` | 219 | 28 | | `james-128-q5` | 307 | 39 |
| `d-james-128-q13` | 285 | 36 | | `james-128-q13` | 337 | 43 |
| `d-james-128-q23` | 335 | 42 | | `james-128-q23` | 353 | 45 |
| `d-james-256-q2` | 390 | 49 | | `james-256-q2` | 578 | 73 |
| `d-james-256-q4` | 446 | 56 | | `james-256-q4` | 618 | 78 |
| `d-james-256-q5` | 479 | 60 | | `james-256-q5` | 637 | 80 |
| `d-james-256-q13` | 633 | 80 | | `james-256-q13` | 711 | 89 |
| `d-james-256-q23` | 738 | 93 | | `james-256-q23` | 756 | 95 |

### 9.3 Public key

A single base-`q` conversion over millions of digits would be quadratic, so
public keys are packed in **groups**: `g` digits into `B` bytes, with `(g, B)`
maximising `g/B` subject to `q^g ≤ 256^B`, searching `B = 1 … 16`.

| q | g | B | efficiency |
|---|---|---|---|
| 2 | 8 | 1 | 100 % |
| 4 | 4 | 1 | 100 % |
| 5 | 31 | 9 | 100.0 % |
| 13 | 28 | 13 | 99.6 % |
| 23 | 7 | 4 | 99.0 % |

The digit stream is, in order:

1. `pk.aa[p][s]` for `p = 0 … n−1`, `s = p … n−1` — that is `n(n+1)/2`
   vectors, at flat index `off[p] + (s − p)` with `off[p] = Σ_{t<p}(n − t)`;
2. `pk.ab[p][s]` for `p = 0 … n−1`, `s = 0 … n_y−1` — `n·n_y` vectors, at flat
   index `p·n_y + s` (D-James only);

each vector contributing its `m` digits in order `0 … m−1`. Full groups are
packed into `B` bytes; the trailing partial group of `rem` digits is packed
into exactly `vec_bytes(rem, q)` bytes, **not** padded to a full group. A
decoder MUST check the total length and reject any group value `≥ q^take`.

### 9.4 Secret key

The 32-byte seed, verbatim.

---

## 10. Test vectors

`djames-py/kat/*.json` pins, per parameter set: the seed, the public key's
length and SHA3-256 digest, and signatures in full for a fixed message list
(`""`, `"abc"`, `bytes(0..31)`, `"D-James known-answer test vector"`). The KAT
seed for a set is `SHAKE256(set name)` truncated to 32 bytes.

These are **regression vectors**: they were produced by the implementations
here and cannot establish agreement with the paper, only agreement between
implementations of this document.

---

## 11. Divergences from ePrint 2026/1650

Three of the paper's statements are self-contradictory. Each is resolved here
in favour of the reading the rest of the paper supports; all three would
benefit from the authors' confirmation.

1. **Dragon index range.** The text names the maps `L_0, …, L_{d−1}` (so
   `k < d`), and §4.2's MinRank argument wants only the first `d` rows of `H̃`
   nonzero — but the displayed sum runs over `q^i ≤ D = q^d + 1`, which
   includes `i = d`. **We use `k ≤ d`**, following the displayed sum and the
   authors' notebook. Note `H̃`'s `λ` block already occupies row `d` via
   `q^d + q^0 ≤ D`, so the "first `d` rows" claim is loose either way.

2. **`q = 4` central monomials.** All four `q = 4` rows (Tables 3, 4, 5, 6)
   give `X³, X¹⁷`. Over `F_4` there is no `X³ = X^(q^i + q^j)`. Every other row
   of every table uses `X² = X^(q⁰+q⁰)` as its low monomial — `F_2` alone uses
   `X³`, because `X²` is `F_2`-linear there and useless as a quadratic term.
   **We read it as `X²`.**

3. **`m` for `q = 4` at 256-bit.** Table 4 lists `m = 180` with `n = 223`,
   `a = 53`, and Table 8 repeats the triple, so it is not an isolated typo. The
   table header defines `m = n − a = 170` and every other row satisfies that
   exactly. **We use `m = 170`**, which resolves the defining relation rather
   than instantiating the literal row; only the authors can say whether `m` or
   `a` is the wrong entry.

Not implemented: `|sig_short|`, the truncated variant, whose verification
requires a hybrid Gröbner solve.

The paper also reports `|sig_fast|` of 334 and 737 bits for `q = 23` where
`⌈n log₂ 23⌉` is 335 and 738; every other entry matches. And its `|PK|` column
reproduces for James but not for D-James — see
[`djames-py/README.md`](djames-py/README.md).

---

## 12. Side channels

Constant-time behaviour is achievable for most of the scheme, and one part of
it is not.

**Public, so unconstrained:** verification in full — the public key, the
signature and the message are all public.

**Constant-time required, and achieved in `djames-rs/`:** base-field and
extension-field arithmetic (no secret-dependent branches, no secret-indexed
table lookups, no division by secrets), the linear algebra over `F_q` that
touches `S`, `T` and `M_Z`, and public-key assembly.

**Not constant-time, inherently:** root finding. The number of roots of the
central polynomial, and the control flow of equal-degree splitting, depend on
secret-derived data. This is a property of HFE inversion, not of this
specification — no HFE implementation known to us avoids it. Implementations
SHOULD keep the arithmetic beneath it constant-time and SHOULD NOT leak
through memory access patterns, but a timing channel on the number of salt
retries and the root count remains.

Rejection sampling leaks only the number of rejections, which is independent
of the accepted values, and the matrix-rejection loops in §6 reject on the
sampled matrix rather than on the accepted one.

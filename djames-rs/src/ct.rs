//! Constant-time primitives.
//!
//! The rule everywhere below the root finder: no branch and no memory index
//! may depend on a secret. These helpers replace the `if` and the `[i]` that
//! would otherwise do so.
//!
//! `wrapping_mul` on small integers is assumed constant time. That holds on
//! every mainstream 64-bit CPU; it is not guaranteed by the language, and on
//! targets with a data-dependent multiplier (some small embedded cores) the
//! carry-less multiply in `gf` would need replacing with a shift-and-mask
//! version. This is the same assumption BearSSL's portable GHASH makes.

/// `0xFF..FF` when `c` is 1, `0` when `c` is 0. Any other input is a bug.
#[inline(always)]
pub fn mask64(c: u64) -> u64 {
    c.wrapping_neg()
}

/// `0xFF..FF` when `a == b`, else 0.
#[inline(always)]
pub fn eq_mask32(a: u32, b: u32) -> u32 {
    let x = a ^ b;
    // x == 0  <=>  (x | -x) has a clear top bit
    let q = (x | x.wrapping_neg()) >> 31;
    (q ^ 1).wrapping_neg()
}

/// `0xFF..FF` when `a != 0`, else 0.
#[inline(always)]
pub fn nz_mask32(a: u32) -> u32 {
    ((a | a.wrapping_neg()) >> 31).wrapping_neg()
}

/// `0xFF..FF` when `a < b` (unsigned), else 0.
#[inline(always)]
pub fn lt_mask32(a: u32, b: u32) -> u32 {
    // Standard borrow-out extraction without a branch.
    let d = (a as u64).wrapping_sub(b as u64);
    ((d >> 63) as u32).wrapping_neg()
}

/// Select `a` if `mask` is all-ones, `b` if it is zero.
#[inline(always)]
pub fn select32(mask: u32, a: u32, b: u32) -> u32 {
    b ^ (mask & (a ^ b))
}

#[inline(always)]
pub fn select64(mask: u64, a: u64, b: u64) -> u64 {
    b ^ (mask & (a ^ b))
}

/// Conditionally swap two slices, in constant time.
#[inline]
pub fn cswap64(mask: u64, a: &mut [u64], b: &mut [u64]) {
    debug_assert_eq!(a.len(), b.len());
    for (x, y) in a.iter_mut().zip(b.iter_mut()) {
        let t = mask & (*x ^ *y);
        *x ^= t;
        *y ^= t;
    }
}

/// `x mod p` for `x < 2^32` and small `p`, without a division instruction.
///
/// Barrett reduction: `mu = floor(2^40 / p)` precomputed by the caller means
/// `x - p * ((x * mu) >> 40)` lands in `[0, 2p)`, and one conditional
/// subtraction finishes it. Integer division by a runtime value is variable
/// time on several architectures, which is why it is avoided here.
#[inline(always)]
pub fn barrett_u32(x: u32, p: u32, mu: u64) -> u32 {
    let q = (((x as u64) * mu) >> 40) as u32;
    let r = x.wrapping_sub(q.wrapping_mul(p));
    // r is now in [0, 2p); subtract p once more if it fits.
    let m = lt_mask32(r, p);
    select32(m, r, r.wrapping_sub(p))
}

/// `floor(2^40 / p)`, the Barrett constant for [`barrett_u32`].
pub const fn barrett_mu(p: u32) -> u64 {
    (1u64 << 40) / (p as u64)
}

/// Overwrite a slice with zeros in a way the optimiser may not remove.
#[inline(never)]
pub fn zeroize64(buf: &mut [u64]) {
    for b in buf.iter_mut() {
        unsafe { core::ptr::write_volatile(b, 0) };
    }
    core::sync::atomic::compiler_fence(core::sync::atomic::Ordering::SeqCst);
}

#[inline(never)]
pub fn zeroize8(buf: &mut [u8]) {
    for b in buf.iter_mut() {
        unsafe { core::ptr::write_volatile(b, 0) };
    }
    core::sync::atomic::compiler_fence(core::sync::atomic::Ordering::SeqCst);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn masks() {
        assert_eq!(eq_mask32(5, 5), u32::MAX);
        assert_eq!(eq_mask32(5, 6), 0);
        assert_eq!(nz_mask32(0), 0);
        assert_eq!(nz_mask32(1), u32::MAX);
        assert_eq!(nz_mask32(u32::MAX), u32::MAX);
        assert_eq!(lt_mask32(3, 4), u32::MAX);
        assert_eq!(lt_mask32(4, 4), 0);
        assert_eq!(lt_mask32(5, 4), 0);
        assert_eq!(lt_mask32(0, u32::MAX), u32::MAX);
        assert_eq!(select32(u32::MAX, 7, 9), 7);
        assert_eq!(select32(0, 7, 9), 9);
    }

    #[test]
    fn barrett_matches_rem() {
        for &p in &[2u32, 3, 4, 5, 13, 23] {
            let mu = barrett_mu(p);
            for x in 0..5000u32 {
                assert_eq!(barrett_u32(x, p, mu), x % p, "x={x} p={p}");
            }
            for x in [65535u32, 100_000, 1_000_000, 16_777_215, 1 << 24] {
                assert_eq!(barrett_u32(x, p, mu), x % p, "x={x} p={p}");
            }
        }
    }

    #[test]
    fn cswap() {
        let (mut a, mut b) = ([1u64, 2, 3], [9u64, 8, 7]);
        cswap64(0, &mut a, &mut b);
        assert_eq!(a, [1, 2, 3]);
        cswap64(u64::MAX, &mut a, &mut b);
        assert_eq!(a, [9, 8, 7]);
        assert_eq!(b, [1, 2, 3]);
    }
}

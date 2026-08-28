//! Keccak-f[1600] and SHAKE256.
//!
//! Written out rather than pulled from a crate so the whole implementation is
//! dependency-free and auditable in one place. The permutation is
//! straight-line: no data-dependent branches, no table lookups indexed by
//! data, so it is constant time with respect to its input.

const RC: [u64; 24] = [
    0x0000_0000_0000_0001,
    0x0000_0000_0000_8082,
    0x8000_0000_0000_808a,
    0x8000_0000_8000_8000,
    0x0000_0000_0000_808b,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8009,
    0x0000_0000_0000_008a,
    0x0000_0000_0000_0088,
    0x0000_0000_8000_8009,
    0x0000_0000_8000_000a,
    0x0000_0000_8000_808b,
    0x8000_0000_0000_008b,
    0x8000_0000_0000_8089,
    0x8000_0000_0000_8003,
    0x8000_0000_0000_8002,
    0x8000_0000_0000_0080,
    0x0000_0000_0000_800a,
    0x8000_0000_8000_000a,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8080,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8008,
];

const RHO: [u32; 24] = [
    1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44,
];

const PI: [usize; 24] = [
    10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1,
];

fn keccak_f1600(a: &mut [u64; 25]) {
    for round in 0..24 {
        // theta
        let mut c = [0u64; 5];
        for x in 0..5 {
            c[x] = a[x] ^ a[x + 5] ^ a[x + 10] ^ a[x + 15] ^ a[x + 20];
        }
        for x in 0..5 {
            let d = c[(x + 4) % 5] ^ c[(x + 1) % 5].rotate_left(1);
            for y in 0..5 {
                a[x + 5 * y] ^= d;
            }
        }
        // rho and pi
        let mut last = a[1];
        for i in 0..24 {
            let j = PI[i];
            let t = a[j];
            a[j] = last.rotate_left(RHO[i]);
            last = t;
        }
        // chi
        for y in 0..5 {
            let row = [
                a[5 * y],
                a[5 * y + 1],
                a[5 * y + 2],
                a[5 * y + 3],
                a[5 * y + 4],
            ];
            for x in 0..5 {
                a[5 * y + x] = row[x] ^ (!row[(x + 1) % 5] & row[(x + 2) % 5]);
            }
        }
        // iota
        a[0] ^= RC[round];
    }
}

const RATE: usize = 136; // SHAKE256: 1600/8 - 2*256/8

/// SHAKE256 as a squeezable sponge.
///
/// `absorb` may be called repeatedly; the first `squeeze` finalises the state,
/// after which further `squeeze` calls continue the same output stream.
#[derive(Clone)]
pub struct Shake256 {
    state: [u64; 25],
    pos: usize,
    squeezing: bool,
}

impl Default for Shake256 {
    fn default() -> Self {
        Self::new()
    }
}

impl Shake256 {
    pub fn new() -> Self {
        Shake256 {
            state: [0u64; 25],
            pos: 0,
            squeezing: false,
        }
    }

    fn xor_byte(&mut self, i: usize, b: u8) {
        self.state[i / 8] ^= (b as u64) << (8 * (i % 8));
    }

    fn byte(&self, i: usize) -> u8 {
        (self.state[i / 8] >> (8 * (i % 8))) as u8
    }

    pub fn absorb(&mut self, data: &[u8]) {
        debug_assert!(!self.squeezing, "absorb after squeeze");
        for &b in data {
            self.xor_byte(self.pos, b);
            self.pos += 1;
            if self.pos == RATE {
                keccak_f1600(&mut self.state);
                self.pos = 0;
            }
        }
    }

    fn finalize(&mut self) {
        self.xor_byte(self.pos, 0x1f); // SHAKE domain separation + first pad bit
        self.xor_byte(RATE - 1, 0x80); // final pad bit
        keccak_f1600(&mut self.state);
        self.pos = 0;
        self.squeezing = true;
    }

    pub fn squeeze(&mut self, out: &mut [u8]) {
        if !self.squeezing {
            self.finalize();
        }
        for b in out.iter_mut() {
            if self.pos == RATE {
                keccak_f1600(&mut self.state);
                self.pos = 0;
            }
            *b = self.byte(self.pos);
            self.pos += 1;
        }
    }
}

/// One-shot SHAKE256.
pub fn shake256(input: &[u8], out: &mut [u8]) {
    let mut s = Shake256::new();
    s.absorb(input);
    s.squeeze(out);
}

/// SHA3-256, used only to digest public keys for the test vectors.
pub fn sha3_256(input: &[u8]) -> [u8; 32] {
    const R: usize = 136; // 1600/8 - 2*256/8
    let mut st = [0u64; 25];
    let mut pos = 0usize;
    let xor = |st: &mut [u64; 25], i: usize, b: u8| {
        st[i / 8] ^= (b as u64) << (8 * (i % 8));
    };
    for &b in input {
        xor(&mut st, pos, b);
        pos += 1;
        if pos == R {
            keccak_f1600(&mut st);
            pos = 0;
        }
    }
    xor(&mut st, pos, 0x06); // SHA-3 domain separation
    xor(&mut st, R - 1, 0x80);
    keccak_f1600(&mut st);
    let mut out = [0u8; 32];
    for (i, o) in out.iter_mut().enumerate() {
        *o = (st[i / 8] >> (8 * (i % 8))) as u8;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(b: &[u8]) -> alloc::string::String {
        use alloc::string::String;
        use core::fmt::Write;
        let mut s = String::new();
        for x in b {
            write!(s, "{x:02x}").unwrap();
        }
        s
    }

    #[test]
    fn shake256_known_answers() {
        // NIST CAVP / standard vectors for the empty input and "abc".
        let mut out = [0u8; 32];
        shake256(b"", &mut out);
        assert_eq!(
            hex(&out),
            "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
        );
        shake256(b"abc", &mut out);
        assert_eq!(
            hex(&out),
            "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739"
        );
    }

    #[test]
    fn squeeze_is_a_stream() {
        // Squeezing in pieces must equal squeezing in one go, across the
        // 136-byte rate boundary.
        let mut a = [0u8; 400];
        shake256(b"stream", &mut a);
        let mut s = Shake256::new();
        s.absorb(b"stream");
        let mut b = [0u8; 400];
        let mut off = 0;
        for chunk in [1usize, 7, 128, 136, 1, 127, 100] {
            let n = chunk.min(400 - off);
            s.squeeze(&mut b[off..off + n]);
            off += n;
        }
        assert_eq!(off, 400);
        assert_eq!(a[..], b[..]);
    }

    #[test]
    fn sha3_256_known_answers() {
        assert_eq!(
            hex(&sha3_256(b"")),
            "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a"
        );
        assert_eq!(
            hex(&sha3_256(b"abc")),
            "3a985da74fe225b2045c172d6bd390bd855f086e3e9d525b46bfe24511431532"
        );
    }
}

type Hash = u64;

/// (2^64-15)/53
/// This has the property that 2^64 mod P = 15, and 15*P < 2^64.
const P: Hash = 348051774975651917;
const R: Hash = 15;
// Largest 64-bit prime
// const P: Hash = 18446744073709551557;
// Largest 64-bit Mersen prime.
// const P: Hash = (1 << 61) - 1;
const BASE: Hash = 256;
pub struct RollingHash {
    sparse_hashes: Vec<Hash>,
}
/// Note that hashes go from back to front to be able to use the low-endian representation.
/// I.e., we compute $\sum_{i=0}^{n-1} BASE^{i} c_{n-1-i} \mod P$.
impl RollingHash {
    /// Stores partial hashes every `s` positions.
    pub fn new(t: &[u8], s: usize) -> Self {
        todo!()
    }
    /// h = (hl * f + hr) % P
    fn mul_add(hl: Hash, f: Hash, hr: Hash) -> Hash {
        let h2 = hl as u128 * f as u128 + hr as u128;
        let high = (h2 >> 64) as u64;
        let low = h2 as u64 % P;
        (high * R + low) % P
    }
    /// Hash `t` one char at a time.
    pub fn linear_simple(t: &[u8]) -> Hash {
        let mut h = 0;
        for &c in t.iter().rev() {
            // h = (h * BASE + c as Hash) % P;
            h = Self::mul_add(h, BASE, c as Hash);
        }
        h
    }
    /// Hash `t` 8 chars at a time.
    // TODO: SIMD-based hashing of 32 chars at a time?
    pub fn linear_fast(t: &[u8]) -> Hash {
        let it = t.chunks_exact(8);
        // First process the offset at the start.
        let mut h: Hash = 0;
        for &c in it.remainder().iter().rev() {
            h = h * BASE + c as Hash;
        }
        // Process remaining chunks of size 8.
        for c in it.rev() {
            // h = (2^64*h + val) % P = R*h+val;
            h *= R;
            let val = u64::from_le_bytes(c.try_into().unwrap());
            let (h2, overflow) = h.overflowing_add(val);
            // In case of overflow, we lost 2^64 which becomes R%P.
            h = h2 + R * overflow as Hash;
            h %= P;
        }
        h
    }
}

#[cfg(test)]
mod test {
    use std::{hint::black_box, num::Wrapping};

    use super::*;
    #[test]
    fn test_remainder() {
        let m = ((1u128 << 64) % P as u128) as Hash;
        assert_eq!(m, R);
    }
    #[test]
    fn random_string() {
        for _ in 0..100 {
            for len in (0..100).chain([100, 1000, 10000].iter().cloned()) {
                let t = (0..len).map(|_| rand::random::<u8>()).collect::<Vec<_>>();
                let h1 = RollingHash::linear_simple(&t);
                let h2 = RollingHash::linear_fast(&t);
                assert_eq!(h1, h2, "Hash mismatch for len = {}", len);
            }
        }
    }
    #[test]
    fn bench_linear_simple() {
        let mut sum = Wrapping(0);
        for _ in 0..100 {
            for len in (0..100).chain([100, 1000, 10000].iter().cloned()) {
                let t = (0..len).map(|_| rand::random::<u8>()).collect::<Vec<_>>();
                sum += RollingHash::linear_simple(&t);
            }
        }
        black_box(sum);
    }
    #[test]
    fn bench_linear_fast() {
        let mut sum = Wrapping(0);
        for _ in 0..100 {
            for len in (0..100).chain([100, 1000, 10000].iter().cloned()) {
                let t = (0..len).map(|_| rand::random::<u8>()).collect::<Vec<_>>();
                sum += RollingHash::linear_fast(&t);
            }
        }
        black_box(sum);
    }
}

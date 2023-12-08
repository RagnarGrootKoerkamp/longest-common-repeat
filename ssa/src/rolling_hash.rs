use std::ops::{Add, Mul, Range, Sub};

/// (2^64-15)/53
/// This has the property that 2^64 mod P = 15, and 15*P < 2^64.
///
/// t.chunks(s).
const P: u64 = 348051774975651917;
const R: u64 = 15;
// Largest 64-bit prime
// const P: Hash = 18446744073709551557;
// Largest 64-bit Mersen prime.
// const P: Hash = (1 << 61) - 1;
const BASE: u64 = 256;

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct Mod(u64);
impl Mul<u64> for Mod {
    type Output = Self;
    fn mul(self, f: u64) -> Self::Output {
        let h2 = f as u128 * self.0 as u128;
        let high = (h2 >> 64) as u64;
        let low = h2 as u64 % P;
        // NOTE: The addition cannot overflow.
        Mod((high * R + low) % P)
    }
}
impl Mul for Mod {
    type Output = Self;
    fn mul(self, f: Self) -> Self::Output {
        self * f.0
    }
}
impl Add for Mod {
    type Output = Self;
    fn add(self, f: Self) -> Self::Output {
        let mut r = self.0 + f.0;
        if r >= P {
            r -= P;
        }
        Mod(r)
    }
}
impl Sub for Mod {
    type Output = Self;
    fn sub(self, f: Self) -> Self::Output {
        self + Mod(P - f.0)
    }
}
impl Mod {
    /// h = (self * f + hr) % P
    fn mul_add(self, f: u64, hr: Mod) -> Self {
        let h2 = self.0 as u128 * f as u128;
        let high = (h2 >> 64) as u64;
        let low = h2 as u64 % P;
        // NOTE: The addition cannot overflow.
        Mod((high * R + low + hr.0) % P)
    }
    /// h = (2^64*h + c) % P = (R*h + c) % P
    fn roll_add(self, c: u64) -> Self {
        let mut h = self.0 * R;
        // In case of overflow, we lost 2^64 which we add back as R.
        let (h2, overflow) = h.overflowing_add(c);
        // NOTE: The addition here cannot overflow.
        h = h2 + R * overflow as u64;
        h %= P;
        Mod(h)
    }
    fn pow(self, mut exp: u64) -> Self {
        if exp == 0 {
            return Mod(1);
        }
        let mut base = self;
        let mut acc = Mod(1);

        while exp > 1 {
            if (exp & 1) == 1 {
                acc = acc * base;
            }
            exp /= 2;
            base = base * base;
        }

        acc * base
    }
}

pub struct RollingHash<'a> {
    text: &'a [u8],
    /// Store a prefix every s positions.
    s: usize,
    log_s: u32,
    /// 256^s
    f: Mod,
    /// 256^-1
    base_inv: Mod,
    // TODO: Also store f^i values?
    prefixes: Vec<Mod>,
}

/// Compute rolling hashes, i.e. $\sum_{i=0}^{n-1} BASE^i c_{i} \mod P$.
///
/// This uses a backwards loop.
impl<'a> RollingHash<'a> {
    /// Stores partial hashes every `s` positions.
    /// Note that `s` must be a power of 2 at least 8 for efficiency.
    /// Answers queries in time `O(min(|t|, |n/s|))`.
    pub fn new(text: &'a [u8], s: usize) -> Self {
        let n = text.len();
        assert!(s % 8 == 0);
        assert!(s.is_power_of_two());
        // p^s
        let f = Mod(R).pow((s / 8) as u64);

        let mut prefixes = Vec::with_capacity(n / s);

        let mut prefix = Mod(0);
        let mut f_acc = Mod(1);
        prefixes.push(prefix);
        for t in text.chunks_exact(s) {
            let h = Self::linear(t);
            prefix = h.mul_add(f_acc.0, prefix);
            prefixes.push(prefix);
            f_acc = f_acc * f;
        }

        Self {
            text,
            s,
            log_s: s.ilog2(),
            f,
            base_inv: Mod(BASE).pow(P - 2),
            prefixes,
        }
    }

    /// Query the hash of the given range.
    /// Compute as
    /// ```txt
    /// |......|..i...|......|...j..|...
    ///        l             r
    /// ======= lookup pl
    /// ===================== lookup pr
    ///        === scan sl   ==== scan sr
    ///           ^ hl           ^ hr
    ///           =============== ans
    ///
    /// ```
    pub fn query(&self, range: Range<usize>) -> Mod {
        let Range { start: i, end: j } = range;
        if range.len() <= 2 * self.s {
            return Self::linear(&self.text[range]);
        }
        let l = i >> self.log_s;
        let r = j >> self.log_s;
        let pl = self.prefixes[l];
        let pr = self.prefixes[r];
        let sl = Self::linear(&self.text[(l << self.log_s)..i]);
        let sr = Self::linear(&self.text[(r << self.log_s)..j]);
        let hl = pl + self.f.pow(l as u64) * sl;
        let hr = pr + self.f.pow(r as u64) * sr;
        (hr - hl) * self.base_inv.pow(i as u64)
    }

    /// Hash `t` 8 chars at a time.
    // TODO: SIMD-based hashing of 32 chars at a time?
    pub fn linear(t: &[u8]) -> Mod {
        let it = t.chunks_exact(8);
        // First process the offset at the start.
        let mut bytes = [0u8; 8];
        bytes[..it.remainder().len()].copy_from_slice(it.remainder());
        let mut h = Mod(u64::from_le_bytes(bytes));

        // Process remaining chunks of size 8.
        for c in it.rev() {
            // h = (2^64*h + val) % P = R*h+val;
            h = h.roll_add(u64::from_le_bytes(c.try_into().unwrap()));
        }
        h
    }

    #[cfg(test)]
    pub fn query_linear(&self, range: Range<usize>) -> Mod {
        Self::linear(&self.text[range])
    }

    /// Hash `t` one char at a time.
    #[cfg(test)]
    fn linear_baseline(t: &[u8]) -> Mod {
        let mut h = Mod(0);
        for &c in t.iter().rev() {
            h = h.mul_add(BASE, Mod(c as u64));
        }
        h
    }
}

#[cfg(test)]
mod test {
    use std::{hint::black_box, num::Wrapping};

    use rand::Rng;

    use super::*;

    #[test]
    fn is_little_endian() {
        #[cfg(not(target_endian = "little"))]
        panic!("Rabin Karp rolling hash implementation only works on little endian machines.");
    }
    #[test]
    fn remainder() {
        let m = ((1u128 << 64) % P as u128) as u64;
        assert_eq!(m, R);
    }
    #[test]
    fn chunked_linear() {
        for _ in 0..100 {
            for len in (0..100).chain([100, 1000, 10000].iter().cloned()) {
                let t = (0..len).map(|_| rand::random::<u8>()).collect::<Vec<_>>();
                let h1 = RollingHash::linear_baseline(&t);
                let h2 = RollingHash::linear(&t);
                assert_eq!(h1, h2, "Hash mismatch for len {len}");
            }
        }
    }
    #[test]
    fn rolling_hash() {
        for s in [8, 16, 32, 64, 128, 1024, 4096] {
            for len in (0..100).chain([100, 1000, 10000].iter().cloned()) {
                let t = (0..len).map(|_| rand::random::<u8>()).collect::<Vec<_>>();
                let rolling_hash = RollingHash::new(&t, s);
                for _ in 0..1000 {
                    let mut i = rand::thread_rng().gen_range(0..=len);
                    let mut j = rand::thread_rng().gen_range(0..=len);
                    if j < i {
                        (i, j) = (j, i);
                    }
                    let h1 = rolling_hash.query(i..j);
                    let h2 = rolling_hash.query_linear(i..j);
                    assert_eq!(h1, h2, "Hash mismatch for len {len}, s {s}, i {i}, j {j}");
                }
            }
        }
    }

    #[test]
    fn bench_linear_simple() {
        let mut sum = Wrapping(0);
        for _ in 0..100 {
            for len in (0..100).chain([100, 1000, 10000].iter().cloned()) {
                let t = (0..len).map(|_| rand::random::<u8>()).collect::<Vec<_>>();
                sum += RollingHash::linear_baseline(&t).0;
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
                sum += RollingHash::linear(&t).0;
            }
        }
        black_box(sum);
    }
}

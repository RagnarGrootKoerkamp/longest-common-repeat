#![allow(unused)]
use std::sync::Mutex;

use assert2::assert;
use fxhash::FxHashSet;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use rdst::RadixSort;

// type B = u8;
// type B = u16;
// type B = u32;
type B = u64;

const H_LOOKUP: [B; 256] = {
    let mut lookup = [1; 256];
    // let a = 0x3c8b_fbb3_95c6_0474u64 as B; // old
    // let a = 0x3c8b_fbb3_95c6_0470u64 as B; // new
    let a = 0x3c8b_fbb3_95c6_0471u64 as B; // newer
    let c = 0x3193_c185_62a0_2b4cu64 as B;
    let g = 0x2032_3ed0_8257_2324u64 as B;
    // let t = 0x2955_49f5_4be2_4456u64 as B; // old
    let t = a ^ c ^ g; // new
    lookup[b'A' as usize] = a;
    lookup[b'C' as usize] = c;
    lookup[b'G' as usize] = g;
    lookup[b'T' as usize] = t;
    lookup[0] = a as B;
    lookup[1] = c as B;
    lookup[2] = g as B;
    lookup[3] = t as B;
    lookup
};

// OLD version
// const DIFF_BITS: u32 = 3;
// const DIFFS: [B; 7] = {
//     let a = H_LOOKUP[b'A' as usize] as B;
//     let c = H_LOOKUP[b'C' as usize] as B;
//     let g = H_LOOKUP[b'G' as usize] as B;
//     let t = H_LOOKUP[b'T' as usize] as B;
//     [0, a ^ c, a ^ g, a ^ t, c ^ g, c ^ t, g ^ t]
// };

// NEW version
const DIFF_BITS: u32 = 2;
const DIFFS: [B; 4] = {
    let a = H_LOOKUP[b'A' as usize] as B;
    let c = H_LOOKUP[b'C' as usize] as B;
    let g = H_LOOKUP[b'G' as usize] as B;
    let t = H_LOOKUP[b'T' as usize] as B;
    [0, a ^ c, a ^ g, a ^ t]
};

const DIFF_CHARS: [(u8, u8); 7] = [
    (b'A', b'A'),
    (b'A', b'C'),
    (b'A', b'G'),
    (b'A', b'T'),
    (b'C', b'G'),
    (b'C', b'T'),
    (b'G', b'T'),
];

fn to_string(bits: usize, s: &mut [u8]) {
    for (i, c) in s.iter_mut().enumerate() {
        *c = b"ACGT"[(bits >> (2 * i)) % 4];
    }
}

fn test_sort(k: u32) -> u64 {
    let mut s = [0; 64];
    let hmax = 1u64 << std::cmp::min(52, B::BITS);
    let mut v: Vec<_> = (0..1u64 << (2 * k))
        .into_par_iter()
        .map(|i| nthash_mask(i as B, k as _))
        .filter(|&h| (h as u64) < hmax)
        .collect();
    eprintln!("sorting len {} ..", v.len());
    v.radix_sort_unstable();
    let mut count = 0u64;
    for [x, y] in v.array_windows() {
        if x == y {
            count += 1;
            // println!("{count:>20}: {x:0width$b}\r", width = B::BITS as usize);
            // panic!();
            // assert!(x != y);
        }
    }
    eprintln!("Count: {count}");
    // assert!(count == 0);
    count
}

fn main() {
    let mut v = vec![];
    // for k in (0..=B::BITS / 2 + 2).step_by(2) {
    for k in 34..=34 {
        eprintln!("k {k} (new)");
        let start = std::time::Instant::now();
        let collisions = search_collision(k);
        v.push((k, collisions.len()));
        eprintln!("Took {:?}", start.elapsed());

        if !collisions.is_empty() {
            let start = std::time::Instant::now();
            resolve_collisions(k, collisions);
            eprintln!("Took {:?}", start.elapsed());
            // return;
        }

        // eprintln!("k {k} (old)");
        // let start = std::time::Instant::now();
        // v.push(test_sort(k));
        // eprintln!("Took {:?}", start.elapsed());
    }
    eprintln!("Counts:");
    for (k, c) in v {
        eprintln!("{k:>2} {c:>20}");
    }
}

/// Returns a list of collisions.
fn search_collision(k: u32) -> Vec<B> {
    let check_bits = std::cmp::min(61, B::BITS);
    let mut max = (1u64 << check_bits).wrapping_sub(1);
    if check_bits == 64 {
        max = u64::MAX;
    }

    let mut collisions = vec![];
    for r in 0..1 << (B::BITS - check_bits) {
        assert!(k % 2 == 0);

        let mut v: Vec<B> = (0..1 << (DIFF_BITS * k / 2))
            .into_par_iter()
            .flat_map_iter(|i| {
                let (hl, hr) = match reduced_hash(k, i) {
                    Some(value) => value,
                    None => return None.into_iter().chain(None),
                };

                let h1 = if ((hl as u64) >> check_bits) == r {
                    Some(hl)
                } else {
                    None
                };
                let h2 = if ((hr as u64) >> check_bits) == r {
                    Some(hr)
                } else {
                    None
                };
                h1.into_iter().chain(h2)
            })
            .collect();
        eprintln!("sorting len {} ..", v.len());
        v.radix_sort_unstable();
        // Iterate over v1 and v2 and find duplicates.
        for &[x, y] in v.array_windows() {
            if x == y && x != 0 {
                eprintln!("Collision at {:0width$b}", x, width = B::BITS as usize);
                collisions.push(x);
            }
        }
    }
    eprintln!("collisions: {}", collisions.len());
    collisions
}

/// Given a collision hash, find the corresponding string.
fn resolve_collisions(k: u32, c: Vec<B>) {
    eprintln!("Resolving collisions:");
    assert!(k % 2 == 0);

    let fragments: Mutex<Vec<(usize, usize)>> = Mutex::new(vec![(usize::MAX, usize::MAX); c.len()]);

    (0..1 << (DIFF_BITS * k / 2)).into_par_iter().for_each(|i| {
        let (hl, hr) = match reduced_hash(k, i) {
            Some(value) => value,
            None => return,
        };

        if c.contains(&hl) {
            let p = c.iter().position(|&x| x == hl).unwrap();
            let mut fragments = fragments.lock().unwrap();
            fragments[p].0 = i;
        }
        if c.contains(&hr) {
            let p = c.iter().position(|&x| x == hr).unwrap();
            let mut fragments = fragments.lock().unwrap();
            fragments[p].1 = i;
        }
    });

    let mut fs = fragments.into_inner().unwrap();

    for (&(il, ir), c) in fs.iter().zip(c) {
        eprintln!(
            "Collision: {:0width$b} {il} {ir}",
            c,
            width = B::BITS as usize
        );
        let cl = get_chars(k, il);
        let cr = get_chars(k, ir);
        let s1 = format!("{}{}", cl.0, cr.0);
        let s2 = format!("{}{}", cl.1, cr.1);

        eprintln!("  s1: {}", s1);
        eprintln!("  s2: {}", s2);
        let s1_hash = nthash_slice(s1.as_bytes());
        let s2_hash = nthash_slice(s2.as_bytes());
        // eprintln!("s1: {s1_hash}");
        // eprintln!("s2: {s2_hash}");
        assert!(s1_hash == s2_hash);
    }
}

fn new_rol(h: B) -> B {
    let h = h.rotate_left(7);
    h ^ (h << 23)
}
fn new_rol_pow(h: B, i: u32) -> B {
    let mut h = h;
    for _ in 0..i {
        h = new_rol(h);
    }
    h
}

fn reduced_hash(k: u32, i: usize) -> Option<(B, B)> {
    let mut h: B = 0;
    for j in 0..k / 2 {
        let c = (i >> (DIFF_BITS * j)) & ((1 << DIFF_BITS) - 1);
        if c >= DIFFS.len() {
            return None;
        }
        h = new_rol(h);
        h ^= unsafe { DIFFS.get_unchecked(c) };
    }
    let hr = h;
    let mut hl = h;
    for j in 0..k / 2 {
        hl = new_rol(hl);
    }
    // eprintln!(
    //     "Hash of {}^{} is {:0width$b}",
    //     get_chars(k, i).0,
    //     get_chars(k, i).1,
    //     h,
    //     width = B::BITS as usize
    // );
    Some((hl, hr))
}

fn get_chars(k: u32, i: usize) -> (String, String) {
    let mut s1 = vec![];
    let mut s2 = vec![];
    for j in 0..k / 2 {
        let c = (i >> (DIFF_BITS * j)) & ((1 << DIFF_BITS) - 1);
        if c >= DIFFS.len() {
            panic!();
        }
        s1.push(DIFF_CHARS[c].0);
        s2.push(DIFF_CHARS[c].1);
    }
    let s1 = String::from_utf8(s1).unwrap();
    let s2 = String::from_utf8(s2).unwrap();
    (s1, s2)
}

#[inline(always)]
fn h_bitpacked(c: B) -> B {
    unsafe { *H_LOOKUP.get_unchecked(c as usize) }
}
#[inline(always)]
fn h_char(c: u8) -> B {
    unsafe { *H_LOOKUP.get_unchecked(c as usize) }
}
#[inline(always)]
pub fn nthash_mask(mut s: B, k: usize) -> B {
    let mut out = 0;
    for idx in 0..k {
        out = new_rol(out);
        out ^= h_bitpacked(s % 4);
        s >>= 2;
    }
    out
}
#[inline(always)]
pub fn nthash_slice(mut s: &[u8]) -> B {
    let k = s.len();
    let mut out = 0;
    for (idx, &c) in s.iter().enumerate() {
        out = new_rol(out);
        out ^= h_char(c);
    }
    out
}

const C: u64 = 0x517cc1b727220a95;
/// Multiply a bitmask by a mixing constant.
pub fn fxhash_mask(s: u64) -> u64 {
    s.wrapping_mul(C)
}
/// Multiply a string of up to 8 chars by a mixing constant.
pub fn fxhash_slice(s: &[u8]) -> u64 {
    let mut arr = [0u8; 8];
    arr[..s.len()].copy_from_slice(s);
    let v = u64::from_be_bytes(arr);
    v.wrapping_mul(C)
}

#[cfg(test)]
mod test {
    use super::*;
    use assert2::assert;

    #[ignore = "changed constants"]
    #[test]
    fn nthash_slice() {
        let s = b"ACGT";
        let h1 = super::nthash_slice(s) as u64;
        let h2 = nthash::ntf64(s, 0, 4);
        if B::BITS == 64 {
            assert!(h1 == h2);
        }
    }

    #[test]
    fn collision_23() {
        let s1 = b"AAGCAACAAAAGAAAGCAAAGAA";
        let s2 = b"CATTCAGAGTCTTTGTGGATTAC";

        let h1 = nthash::ntf64(s1, 0, 23);
        let h2 = nthash::ntf64(s2, 0, 23);
        assert!(h1 == h2);
    }

    #[test]
    fn collision_24() {
        let s1 = b"CAAGAAAGAAACACCAACAAACAG";
        let s2 = b"GCCTCCGTAAGTCTGTCGTCAGAT";

        let h1 = nthash::ntf64(s1, 0, 24);
        let h2 = nthash::ntf64(s2, 0, 24);
        assert!(h1 == h2);
    }

    #[test]
    fn test_invertible() {
        let a = h_char(b'A');
        let c = h_char(b'C');
        let g = h_char(b'G');
        let t = h_char(b'T');
        assert!(a ^ c ^ g ^ t == 0);
        let mut v = [0; B::BITS as usize];
        let mut x = a ^ c;
        let mut y = a ^ g;
        for i in 0..B::BITS / 2 {
            v[2 * i as usize] = x;
            v[2 * i as usize + 1] = y;
            x = new_rol(x);
            y = new_rol(y);
        }
        assert!(gaussian_elimination(v));
    }

    fn gaussian_elimination(mut words: [B; B::BITS as usize]) -> bool {
        let mut inv = true;
        // i: fixed rows.
        let mut i = 0;
        // j: fixed columns.
        for j in 0..B::BITS as usize {
            // Find a pivot row k with bit i set.
            let mut k = i;
            while k < B::BITS as usize && words[k] & (1 << j) == 0 {
                k += 1;
            }
            if k == B::BITS as usize {
                inv = false;
                continue;
            }
            // Make row i the pivot row.
            words.swap(i, k);
            for k in (i + 1)..B::BITS as usize {
                if words[k] & (1 << j) != 0 {
                    words[k] ^= words[i];
                }
            }
            i += 1;
        }
        for w in words {
            eprintln!("{:0width$b}", w, width = B::BITS as usize);
        }
        inv
    }

    #[test]
    fn test_leading_zeros() {
        // let a = h_char(b'A');
        // let c = h_char(b'C');
        // let g = h_char(b'G');
        // let t = h_char(b'T');
        let a = rand::random::<u64>();
        let c = rand::random::<u64>();
        let g = rand::random::<u64>();
        let t = a ^ c ^ g;
        assert!(a ^ c ^ g ^ t == 0);
        let chars = [a, c, g, t];
        for k in (1..64).step_by(2) {
            for c1 in 0..4 {
                for c2 in 0..4 {
                    let h = chars[c1] ^ chars[c2].rotate_left(7 * k);
                    // let h = chars[c1] ^ new_rol_pow(chars[c2], k);
                    let lz = h.leading_zeros();
                    if lz > 5 {
                        eprintln!("k={k:2} c1={c1} c2={c2} h={h:064b} lz={lz}");
                    }
                }
            }
        }
    }
}

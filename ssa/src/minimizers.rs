use std::collections::VecDeque;

use crate::rolling_hash::RollingHash;

/// Finds the position of the length-k minimizer in each window of length w.
/// NOTE: w is the number of kmers from which we find a minimizer, in a span of k+w-1 bases total.
pub fn minimizers(t: &[u8], k: usize, w: usize) -> Vec<usize> {
    assert!(w >= 1);
    // 1. Make kmer windows.
    // 2. Compute (rolling) hash of each kmer.
    // 3. Make w windows.
    // 4. Find min hash in each w window.

    let mut minimizers = Vec::new();
    let mut stack = VecDeque::new();

    for (i, h) in RollingHash::kmer_hashes_rev(t, k) {
        // 1. drop worse minimizers.
        while let Some((_, h2)) = stack.back() {
            if h < *h2 {
                stack.pop_back();
            } else {
                break;
            }
        }
        // 2. push new minimizer.
        stack.push_back((i, h));

        // 3. Mark front as best, if not pushed already.
        let i2 = stack.front().unwrap().0;
        assert!(i2 <= i + w);
        if i + w <= t.len() && minimizers.last() != Some(&i2) {
            minimizers.push(i2);
        }
        // 4. Pop the front if it moves out of the window.
        if i2 == i + w - 1 {
            stack.pop_front();
        }
    }
    minimizers.reverse();
    minimizers
}

/// Finds the position of the length-k minimizer in each window of length w.
/// NOTE: w is the number of kmers from which we find a minimizer, in a span of k+w-1 bases total.
pub fn minimizers_nthash(t: &[u8], k: usize, w: usize) -> Vec<usize> {
    assert!(w >= 1);
    // 1. Make kmer windows.
    // 2. Compute (rolling) hash of each kmer.
    // 3. Make w windows.
    // 4. Find min hash in each w window.

    let mut minimizers = Vec::new();
    let mut stack = VecDeque::new();

    for (i, h) in nthash::NtHashForwardIterator::new(t, k)
        .unwrap()
        .enumerate()
    {
        // for (i, h) in RollingHash::kmer_hashes_rev(t, k) {
        // 1. drop worse minimizers.
        while let Some((_, h2)) = stack.back() {
            if h <= *h2 {
                stack.pop_back();
            } else {
                break;
            }
        }
        // 2. push new minimizer.
        stack.push_back((i, h));

        // 3. Mark front as best, if not pushed already.
        let i2 = stack.front().unwrap().0;
        assert!(i <= i2 + w);
        if i >= w - 1 && minimizers.last() != Some(&i2) {
            minimizers.push(i2);
        }
        // 4. Pop the front if it moves out of the window.
        if i == i2 + w - 1 {
            stack.pop_front();
        }
    }
    minimizers
}

/// ntHash constants.
static LUT: [u64; 128] = {
    let mut l = [0u64; 128];
    l[b'A' as usize] = 0x3c8bfbb395c60474u64;
    l[b'C' as usize] = 0x3193c18562a02b4cu64;
    l[b'G' as usize] = 0x20323ed082572324u64;
    l[b'T' as usize] = 0x295549f54be24456u64;
    l
};

const C: u64 = 0x517cc1b727220a95;
fn get(c: u8) -> u64 {
    // c as u64 * C
    unsafe { *LUT.get_unchecked(c as usize) }
}

pub fn minimizers_daniel(s: &[u8], k: usize, w: usize) -> Vec<usize> {
    let w = w + k - 1;
    let mut v = vec![];
    minimizers_daniel_inner(s, w, k, 0, |i, m| {
        v.push(i);
    });
    v
}

/// Robust winnowing.
/// https://gist.github.com/Daniel-Liu-c0deb0t/7078ebca04569068f15507aa856be6e8
/// Fixed 2 bugs: h<=min instead of h<min, and i+1 instead of i.
#[inline(always)]
fn minimizers_daniel_inner(s: &[u8], w: usize, k: usize, prev: u64, mut f: impl FnMut(usize, u64)) {
    assert!(k <= w);
    let mut min = 0;
    let mut min_idx = 0;
    let mut curr = 0;

    for (i, win) in s.windows(w).enumerate() {
        if i == 0 || i > min_idx {
            let (m_idx, m, c) = minimum(win, k, prev, i);
            min_idx = i + m_idx;
            min = m;
            curr = c;
            f(min_idx, min);
        } else {
            curr =
                curr.rotate_left(1) ^ get(win[w - 1 - k]).rotate_left(k as u32) ^ get(win[w - 1]);
            let h = prev ^ curr;

            if h <= min {
                min_idx = i + w - k;
                min = h;
                f(min_idx, min);
            }
        }
    }
}

#[inline(always)]
/// Get the rightmost minimum kmer.
fn minimum(s: &[u8], k: usize, prev: u64, offset: usize) -> (usize, u64, u64) {
    let mut curr = 0;

    for (i, &b) in s[..k].iter().enumerate() {
        curr ^= get(b).rotate_left((k - 1 - i) as u32);
    }

    let mut min = prev ^ curr;
    let mut min_idx = 0;

    for (i, &b) in s[k..].iter().enumerate() {
        curr = curr.rotate_left(1) ^ get(s[i]).rotate_left(k as u32) ^ get(b);
        let h = prev ^ curr;

        if h <= min {
            min = h;
            min_idx = i + 1;
        }
    }

    (min_idx, min, curr)
}

#[cfg(test)]
mod test {
    #[test]
    fn small() {
        let t = &b"LRagnarWasHerezWI8sOLRagnarWasHereLr"[..];
        let o0 = 1;
        let o1 = 21;
        let secret_len = 13;
        for w in 1..13 {
            for k in 1..13 {
                if w + k > secret_len - 2 {
                    continue;
                }
                let ms = super::minimizers(t, k, w);
                let mut ok = false;
                for &m in &ms {
                    if o0 <= m && m <= o0 + secret_len - k && ms.contains(&(m + o1 - o0)) {
                        ok = true;
                        break;
                    }
                }
                assert!(
                    ok,
                    "did not find two corresponding minimizers for k={k} w={w}\n{ms:?}"
                );
            }
        }
    }

    fn verify_minimizers(w: usize, ms: &Vec<usize>) {
        if ms.is_empty() {
            return;
        }
        for i in 0..ms.len() - 1 {
            // Increasing
            assert!(
                ms[i] <= ms[i + 1],
                "Error at position {i}: {} <= {}",
                ms[i],
                ms[i + 1]
            );
            // Difference at most w.
            assert!(
                ms[i + 1] - ms[i] <= w,
                "Gap too large {} {} w={w}",
                ms[i],
                ms[i + 1]
            );
        }
    }

    fn test_minimizers(f: impl Fn(&[u8], usize, usize) -> Vec<usize>) {
        for n in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000, 10000] {
            for w in 1..10 {
                for k in 1..=2 * w {
                    if k > n {
                        continue;
                    }
                    let t = (0..n)
                        .map(|_| b"ACTG"[rand::random::<usize>() % 4])
                        .collect::<Vec<_>>();
                    let ms = f(&t, k, w);
                    verify_minimizers(w, &ms);
                }
            }
        }
    }

    #[test]
    fn minimizers_naive() {
        test_minimizers(super::minimizers);
    }
    #[test]
    fn minimizers_daniel() {
        test_minimizers(super::minimizers_daniel);
    }
    #[test]
    fn minimizers_nthash() {
        test_minimizers(super::minimizers_nthash);
    }
    #[test]
    fn minimizers_nthash_same() {
        test_minimizers(|t, k, w| {
            let ms1 = super::minimizers_nthash(t, k, w);
            let ms2 = super::minimizers_daniel(t, k, w);
            assert_eq!(
                ms1,
                ms2,
                "Difference for {} {} {}\n{}",
                t.len(),
                k,
                w,
                std::str::from_utf8(t).unwrap(),
            );
            ms1
        });
    }
}

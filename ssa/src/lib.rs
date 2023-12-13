#![feature(impl_trait_in_assoc_type, slice_group_by)]

pub mod rolling_hash;

use rdst::{RadixKey, RadixSort};
use rolling_hash::Mod;

use crate::rolling_hash::RollingHash;

/// A sparse suffix array.
/// Implementation based on https://arxiv.org/pdf/2310.09023.pdf
/// C++ code at https://github.com/lorrainea/SSA/blob/main/PA/ssa.cc
pub struct Ssa {
    pub sa: Vec<usize>,
    pub lcp: Vec<usize>,
}

#[derive(Debug, Clone, Copy)]
struct IH {
    /// Index in the original text.
    idx: usize,
    /// Hash of a prefix-extension of the suffix starting at `i`.
    h: Mod,
}

/// Sort `IH` by the contained hash.
impl RadixKey for IH {
    const LEVELS: usize = 8;
    #[inline]
    fn get_level(&self, level: usize) -> u8 {
        (self.h.0 >> (level * 8)) as u8
    }
}

fn group_len(v: &mut [IH], i: usize) -> usize {
    v[i..].iter().take_while(|x| x.h == v[i].h).count()
}

impl Ssa {
    pub fn new(t: &[u8], idxs: Vec<usize>) -> Self {
        assert!(!idxs.is_empty());
        let n = t.len();
        let max_l = 1 << n.ilog2();
        let b = idxs.len();
        let s = if n == 1 {
            8
        } else {
            (n / n.ilog2() as usize).next_power_of_two().max(8)
        };
        let hasher = rolling_hash::RollingHash::new(t, s);

        let mut starts = idxs
            .into_iter()
            .map(|idx| IH {
                idx,
                h: Mod::default(),
            })
            .collect::<Vec<_>>();

        let mut lcp = vec![usize::MAX; b - 1];

        // Contains [witness, length, idxs...] for each group consecutively.
        let mut cache = vec![];

        fn witness(idx: usize, t: &[u8], cache: &Vec<usize>) -> usize {
            if idx < t.len() {
                idx
            } else {
                cache[idx - t.len()]
            }
        }

        /// Given a slice of indices that already have the given lcp, sort them in-place and write the LCP array.
        fn dfs(
            l: usize,
            group_lcp: usize,
            t: &[u8],
            hasher: &RollingHash,
            starts: &mut [IH],
            cache: &mut Vec<usize>,
            lcp_out: &mut [usize],
        ) {
            let n = starts.len();
            assert!(lcp_out.len() == n - 1);
            if n <= 1 {
                return;
            }

            if l == 0 {
                // Simply sort the groups by their next character.
                lcp_out.fill(group_lcp);

                for IH { idx, h } in starts.iter_mut() {
                    let idx = witness(*idx, t, cache);
                    *h = hasher.query(idx + group_lcp..idx + group_lcp + 1);
                }
                starts.radix_sort_unstable();
                return;
            }

            assert!(l.is_power_of_two(), "{l} is not a power of two");

            // TODO: Resolve size-2 groups using LCP.

            // Group start..end by their 2^len prefix extension.
            // First, compute hashes.
            for IH { idx, h } in starts.iter_mut() {
                let idx = witness(*idx, t, cache);
                *h = hasher.query(idx + group_lcp..idx + group_lcp + l);
            }
            // Second, sort by hashes.
            starts.radix_sort_unstable();
            // Third, count groups.
            let num_groups = starts.group_by(|a, b| a.h == b.h).count();
            // Fourth, recurse into groups.
            if num_groups == 1 {
                // One big group: Recurse with increased LCP length.
                dfs(l / 2, group_lcp + l, t, hasher, starts, cache, lcp_out);
                return;
            }
            if num_groups == n {
                // All groups are singletons: Recurse with original LCP length.
                dfs(l / 2, group_lcp, t, hasher, starts, cache, lcp_out);
                return;
            }
            // Otherwise:
            // 1. Move groups > 1 to the cache.
            // 2. Replace them by witness values in the main array.
            // 3. Recursively sort the main array.
            // 4. Insert groups back into the main array.
            // 5. Recurse on the groups.

            let old_cache_len = cache.len();

            // 1
            let mut i = 0;
            let mut j = 0;
            while i < n {
                let group_len = group_len(starts, i);
                assert!(group_len > 0);
                if group_len == 1 {
                    starts[j] = starts[i];
                } else {
                    let group_idx = t.len() + cache.len();
                    let witness = witness(starts[i].idx, t, cache);
                    // 1. Move groups to the cache.
                    cache.push(witness);
                    cache.push(group_len);
                    cache.extend(starts[i..i + group_len].iter().map(|x| x.idx));
                    // 2. Write sentinel to main array.
                    starts[j].idx = group_idx;
                }
                i += group_len;
                j += 1;
            }
            assert_eq!(i, n);
            assert_eq!(j, num_groups);

            // 3. Recursively sort the main array.
            dfs(
                l / 2,
                group_lcp,
                t,
                hasher,
                &mut starts[..j],
                cache,
                &mut lcp_out[..j - 1],
            );

            // 4. Insert cached groups back into the main array.
            while j > 0 {
                j -= 1;
                if starts[j].idx < t.len() + old_cache_len {
                    i -= 1;
                    starts[i].idx = starts[j].idx;
                    starts[i].h = Mod::NONE;
                } else {
                    let cache_idx = starts[j].idx - t.len();
                    let group_len = cache[cache_idx + 1];
                    i -= group_len;
                    for k in 0..group_len {
                        starts[i + k].idx = cache[cache_idx + 2 + k];
                        starts[i + k].h = Mod(cache_idx as u64);
                    }
                }
                if i > 0 {
                    lcp_out[i - 1] = lcp_out[j - 1];
                }
            }
            assert_eq!(i, 0);

            cache.shrink_to(old_cache_len);

            // 5. Recurse on the groups.
            while i < n {
                if starts[i].h == Mod::NONE {
                    i += 1;
                    continue;
                }
                let group_len = group_len(starts, i);
                dfs(
                    l / 2,
                    group_lcp + l,
                    t,
                    hasher,
                    &mut starts[i..i + group_len],
                    cache,
                    &mut lcp_out[i..i + group_len - 1],
                );
                i += group_len;
            }
        }

        dfs(max_l, 0, t, &hasher, &mut starts, &mut cache, &mut lcp);

        Self {
            sa: starts.iter().map(|x| x.idx).collect(),
            lcp,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn lcp(t: &[u8], a: usize, b: usize) -> usize {
        std::iter::zip(&t[a..], &t[b..])
            .take_while(|(a, b)| a == b)
            .count()
    }

    fn verify(ssa: Ssa, t: &[u8]) {
        let b = ssa.sa.len();
        assert_eq!(ssa.lcp.len(), b - 1);
        for i in 0..b - 1 {
            assert!(
                t[ssa.sa[i]..] <= t[ssa.sa[i + 1]..],
                "Bad order at position {i}"
            );
            assert_eq!(
                lcp(&t, ssa.sa[i], ssa.sa[i + 1]),
                ssa.lcp[i],
                "Bad LCP at position {i}. Got {}",
                ssa.lcp[i]
            );
        }
    }

    #[test]
    fn small() {
        for t in [
            // &b"aaaa"[..],
            // &b"aa"[..],
            // &b"baa"[..],
            // &b"abracadabra"[..],
            // &[4, 2, 1, 1, 2],
            &[0, 0, 0, 0],
        ] {
            let idxs = (0..t.len()).collect::<Vec<_>>();
            let ssa = Ssa::new(t, idxs);
            verify(ssa, t);
        }
    }

    #[test]
    fn random() {
        for &fraction in [0.01, 0.1, 0.5, 0.9, 0.99, 1.0].iter().rev() {
            for len in (1..100).chain([100, 1000, 10000, 100000].iter().copied()) {
                let t = (0..len)
                    .map(|_| rand::random::<u8>() % 4)
                    .collect::<Vec<_>>();
                let idxs = (0..len)
                    .filter(|_| rand::random::<f64>() < fraction)
                    .collect::<Vec<_>>();
                if idxs.is_empty() {
                    continue;
                }
                let ssa = Ssa::new(&t, idxs);
                verify(ssa, &t);
            }
        }
    }
}

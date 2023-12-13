#![feature(slice_partition_dedup)]
use std::time::Duration;

use rdst::RadixSort;

// Group using a hashmap.
fn group_hash(v: &[u64]) -> (usize, Duration) {
    let start = std::time::Instant::now();
    let mut map = std::collections::HashMap::new();
    for &x in v {
        *map.entry(x).or_insert(0) += 1;
    }
    (map.len(), start.elapsed())
}

// Group using radix sort.
fn group_radix(v: &mut [u64]) -> (usize, Duration) {
    let start = std::time::Instant::now();
    v.radix_sort_unstable();
    let len = v.partition_dedup().0.len();
    (len, start.elapsed())
}

fn main() {
    let total = 100000000;
    for n in [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000] {
        let mut t = (0..total)
            .map(|_| rand::random::<u64>() % (n as u64 / 2))
            .collect::<Vec<_>>();

        let mut ht = Duration::default();
        let mut rt = Duration::default();

        for start in (0..total).step_by(n) {
            let h = group_hash(&t[start..start + n]);
            let r = group_radix(&mut t[start..start + n]);
            ht += h.1;
            rt += r.1;
            assert_eq!(h.0, r.0);
        }
        eprintln!(
            "n = {}\n  hash: {:>8?}\n radix: {:>8?}",
            n,
            ht / (total / n) as u32,
            rt / (total / n) as u32
        );
    }
}

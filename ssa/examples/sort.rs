#![feature(slice_partition_dedup)]
use std::time::Duration;

use rdst::RadixSort;

// Group using a hashmap.
fn group_hash(v: &[u64]) -> (usize, Duration) {
    let start = std::time::Instant::now();
    let mut map = std::collections::HashMap::with_capacity(v.len());
    for &x in v {
        *map.entry(x).or_insert(0) += 1;
    }
    (map.len(), start.elapsed())
}

// Group using radix sort.
fn group_radix(v: &mut [u64]) -> (usize, Duration) {
    let start = std::time::Instant::now();
    v.radix_sort_builder().with_single_threaded_tuner().sort();
    let len = v.partition_dedup().0.len();
    (len, start.elapsed())
}

// Group using standard sort.
fn group_sort(v: &mut [u64]) -> (usize, Duration) {
    let start = std::time::Instant::now();
    v.sort();
    let len = v.partition_dedup().0.len();
    (len, start.elapsed())
}

fn main() {
    let total = 100000000;
    for n in [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000] {
        let t = (0..total)
            // .map(|_| rand::random::<u64>() % (n as u64 / 2))
            .map(|_| rand::random::<u64>() % (n as u64 * 2))
            // .map(|_| rand::random::<u64>())
            .collect::<Vec<_>>();

        let mut ht = Duration::default();
        let mut rt = Duration::default();
        let mut st = Duration::default();

        for start in (0..total).step_by(n) {
            let v = &mut t[start..start + n].to_vec();
            let h = group_hash(v);
            let v = &mut t[start..start + n].to_vec();
            let r = group_radix(v);
            let v = &mut t[start..start + n].to_vec();
            let s = group_sort(v);
            ht += h.1;
            rt += r.1;
            st += s.1;
            assert_eq!(h.0, r.0);
            assert_eq!(h.0, s.0);
        }
        eprintln!(
            "n = {}\n  hash: {:>8?}\n radix: {:>8?}\n  sort: {:>8?}",
            n,
            ht / (total / n) as u32,
            rt / (total / n) as u32,
            st / (total / n) as u32
        );
    }
}

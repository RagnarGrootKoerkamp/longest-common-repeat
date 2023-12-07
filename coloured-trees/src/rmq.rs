use std::{default::Default, ops::Range};

/// Range Minimum Query.
pub trait Rmq<T: Ord> {
    fn new(a: &[T]) -> Self;
    fn query(&self, range: Range<usize>) -> T;
}

/// O(1) query, O(n lg n) words of space.
pub struct SparseTable<T> {
    /// The number of elements.
    n: usize,
    /// Packed rows of n elements.
    table: Vec<T>,
}

impl<T: Ord + Clone + Copy + Default> Rmq<T> for SparseTable<T> {
    fn new(a: &[T]) -> Self {
        let n = a.len();
        let logn = n.ilog2() as usize;
        let mut table = vec![T::default(); n * (logn + 1)];
        table[..n].copy_from_slice(a);
        for k in 1..=logn {
            let len = 1 << (k - 1);
            for i in 0..n - len {
                table[k * n + i] = T::min(table[(k - 1) * n + i], table[(k - 1) * n + i + len]);
            }
        }

        Self { n, table }
    }

    /// Query the minimum of [l, r).
    /// 0-based, right exclusive.
    fn query(&self, range: Range<usize>) -> T {
        let Range { start: l, end: r } = range;
        let k = (r - l).ilog2() as usize;
        T::min(
            self.table[k * self.n + l],
            self.table[k * self.n + r - (1 << k)],
        )
    }
}

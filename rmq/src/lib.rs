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

/// O(1) query, O(n) words of space.
/// Taken from  https://codeforces.com/blog/entry/78931
pub struct MaskRmq<T> {
    /// The input values.
    a: Vec<T>,
    /// A sparse table on the minima of blocks of size W.
    sparse: SparseTable<T>,
    /// For a position i, consider a[i..i+W].
    /// The mask at i has bit j set if a[i+j']>a[i+j] for all j'<j.
    masks: Vec<u64>,
}

type B = u64;
const W: usize = B::BITS as usize;

impl<T: Ord + Clone + Copy + Default> Rmq<T> for MaskRmq<T> {
    fn new(a: &[T]) -> Self {
        let chunk_mins: Vec<_> = a
            .chunks(W)
            .map(|block| *block.iter().min().unwrap())
            .collect();
        let sparse = SparseTable::new(&chunk_mins);

        let mut masks = vec![0; a.len()];
        let mut mask: B = 0;
        for (i, &x) in a.iter().enumerate().rev() {
            mask <<= 1;
            while mask > 0 {
                // Clear bits of values larger than the current one.
                if a[i + mask.trailing_zeros() as usize] > x {
                    // Clear the lsb.
                    mask &= mask - 1;
                } else {
                    break;
                }
            }
            mask |= 1;
            masks[i] = mask;
        }
        Self {
            a: a.to_vec(),
            sparse,
            masks,
        }
    }

    fn query(&self, range: Range<usize>) -> T {
        let mask_min = |mask: B| -> T {
            let offset = W - 1 - mask.leading_zeros() as usize;
            self.a[range.start + offset]
        };
        if range.len() <= W {
            return mask_min(self.masks[range.start] & ((1 << range.len()) - 1));
        }
        let ends = T::min(
            // head
            mask_min(self.masks[range.start]),
            // tail
            mask_min(self.masks[range.end - W]),
        );

        let blocks = self
            .sparse
            .query(range.start.div_ceil(W)..range.end.div_ceil(W));

        T::min(ends, blocks)
    }
}

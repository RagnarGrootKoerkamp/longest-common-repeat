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
        assert!(!range.is_empty());
        let mask_min = |pos: usize, bitmask: u64| -> T {
            let mask = self.masks[pos] & bitmask;
            let offset = W - 1 - mask.leading_zeros() as usize;
            self.a[pos + offset]
        };
        if range.len() < W {
            let bitmask = (1u64 << (range.len())).wrapping_sub(1);
            return mask_min(range.start, bitmask);
        }
        let ends = T::min(
            // head
            mask_min(range.start, B::MAX),
            // tail
            mask_min(range.end - W, B::MAX),
        );

        let range = range.start.div_ceil(W)..range.end / W;
        if range.is_empty() {
            return ends;
        }
        let blocks = self.sparse.query(range);
        T::min(ends, blocks)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn test_rmq<R: Rmq<u32>>() {
        for n in (1..10).chain((100..1000).step_by(200)) {
            let a = (0..n).map(|_| rand::random::<u32>()).collect::<Vec<_>>();
            let rmq = R::new(&a);
            for i in 0..n {
                for j in i + 1..=n {
                    let test_ans = rmq.query(i..j);
                    let real_ans = *a[i..j].iter().min().unwrap();
                    assert_eq!(
                        test_ans,
                        real_ans,
                        "Failure for n={n} i={i} j={j} len={}.",
                        j - i,
                    );
                }
            }
        }
    }

    #[test]
    fn sparse_table() {
        test_rmq::<SparseTable<u32>>();
    }
    #[test]
    fn mask_rmq() {
        test_rmq::<MaskRmq<u32>>();
    }
}

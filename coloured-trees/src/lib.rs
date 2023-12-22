use itertools::Itertools;
use rmq::Rmq;
use std::{cmp::max, collections::HashMap};

type Node = usize;
type Weight = usize;

#[derive(Debug)]
pub struct Tree {
    /// Leaf labels.
    pub sa: Vec<Node>,
    /// LCPs between adjacent leafs.
    /// NOTE: This must have the same length as leafs and be end-padded with a 0.
    pub lcp: Vec<Weight>,
}

/// Find a pair of nodes (u, v) such that LCP_a(u, v) + LCP_b(u, v) is maximized.
/// The trees must have the same set of nodes.
pub fn max_common_weight(a: &Tree, b: &Tree) -> (Weight, (Node, Node)) {
    assert_eq!(a.lcp.len(), a.sa.len());
    assert_eq!(a.lcp.last(), Some(&0));
    assert_eq!(b.lcp.len(), b.sa.len());
    assert_eq!(b.lcp.last(), Some(&0));
    assert_eq!(a.sa.len(), b.sa.len());

    let mut b_idx = get_permutation(a, b);

    let mut ans = (0, (Node::MAX, Node::MAX));

    let rmq = rmq::MaskRmq::new(&b.lcp);

    // Inclusive start pos in a of range of subtree, and right-lcp.
    let mut stack = vec![(0, 0)];
    for (i, &a_lcp_right) in a.lcp.iter().enumerate() {
        let mut start = i;
        while !stack.is_empty() && stack.last().unwrap().1 >= a_lcp_right {
            let (old_start, a_lcp) = stack.pop().unwrap();

            // Merge ranges.
            // TODO: use a faster algorithm. This is only nice if trees have depth O(log n).
            b_idx[old_start..=i].sort();

            // Update answer for all adjacent pairs from opposite ranges.
            // TODO: When ranges are very unbalanced, only check the insertion points.
            for (&(bl, al), &(br, ar)) in b_idx[old_start..=i].iter().tuple_windows() {
                if (al < start) ^ (ar < start) {
                    let b_lcp = rmq.query(bl..br);
                    ans = max(ans, (a_lcp + b_lcp, (al, ar)));
                }
            }

            start = old_start;
        }
        stack.push((start, a_lcp_right));
    }

    ans
}

/// Find the permutation from a nodes to b nodes.
/// Returns a vec of (b_idx, a_idx) pairs.
fn get_permutation(a: &Tree, b: &Tree) -> Vec<(usize, usize)> {
    let mut b_inv: HashMap<usize, usize> = HashMap::new();
    for (i, &node) in b.sa.iter().enumerate() {
        b_inv.insert(node, i);
    }
    a.sa.iter()
        .enumerate()
        .map(|(i, ai)| (b_inv[&ai], i))
        .collect()
}

#[cfg(test)]
fn max_common_weight_naive(a: &Tree, b: &Tree) -> (Weight, (Node, Node)) {
    use std::cmp::min;

    let mut ans = (0, (Node::MAX, Node::MAX));
    let a_rmq = rmq::MaskRmq::new(&a.lcp);
    let b_rmq = rmq::MaskRmq::new(&b.lcp);
    let p = get_permutation(a, b);
    for i in 0..a.sa.len() {
        for j in i + 1..a.sa.len() {
            let u = p[i].0;
            let v = p[j].0;
            let lcp = a_rmq.query(i..j);
            let lcp2 = b_rmq.query(min(u, v)..max(u, v));
            ans = max(ans, (lcp + lcp2, (i, j)));
        }
    }
    ans
}

#[cfg(test)]
mod test {
    use rand::seq::SliceRandom;

    use super::*;
    #[test]
    fn max_common_weight_small() {
        let a = Tree {
            sa: vec![0, 1, 2, 3, 4, 5, 6, 7],
            lcp: vec![2, 3, 1, 4, 2, 1, 3, 0],
        };
        let b = Tree {
            sa: vec![7, 6, 5, 4, 3, 2, 1, 0],
            lcp: vec![0, 1, 3, 1, 1, 2, 1, 0],
        };
        max_common_weight(&a, &b);
    }

    #[test]
    fn random() {
        for n in 1..300 {
            let mut a = Tree {
                sa: (0..n).collect(),
                lcp: (0..n).map(|_| rand::random::<usize>() % 10).collect(),
            };
            a.lcp[n - 1] = 0;
            let mut b_sa = (0..n).collect_vec();
            b_sa.shuffle(&mut rand::thread_rng());
            let mut b = Tree {
                sa: b_sa,
                lcp: (0..n).map(|_| rand::random::<usize>() % 10).collect(),
            };
            b.lcp[n - 1] = 0;

            let ans = max_common_weight(&a, &b);
            let naive = max_common_weight_naive(&a, &b);
            assert_eq!(
                ans.0, naive.0,
                "Failure at n={n}\n{a:?}\n{b:?}\nans  : {ans:?}\nnaive: {naive:?}",
            );
        }
    }
}

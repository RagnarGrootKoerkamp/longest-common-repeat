use itertools::Itertools;
use rmq::Rmq;
use std::cmp::max;

type Node = usize;
type Weight = u64;

#[derive(Debug)]
pub struct Tree {
    /// Permutation of the leafs. Must be 0..n-1
    leafs: Vec<Node>,
    /// LCPs between adjacent leafs.
    /// NOTE: This must have the same length as leafs and be end-padded with a 0.
    lcp: Vec<Weight>,
}

/// Find a pair of nodes (u, v) such that LCP_a(u, v) + LCP_b(u, v) is maximized.
pub fn max_common_weight(a: &Tree, b: &Tree) -> (Weight, (Node, Node)) {
    assert_eq!(a.lcp.len(), a.leafs.len());
    assert_eq!(a.lcp.last(), Some(&0));
    assert_eq!(b.lcp.len(), b.leafs.len());
    assert_eq!(b.lcp.last(), Some(&0));
    assert_eq!(a.leafs.len(), b.leafs.len());

    let mut b_inv = vec![0; b.leafs.len()];
    for (i, &node) in b.leafs.iter().enumerate() {
        b_inv[node] = i;
    }
    let mut b_idx: Vec<_> = a.leafs.iter().map(|&i| (b_inv[i], i)).collect();
    drop(b_inv);

    let mut ans = (0, (Node::MAX, Node::MAX));

    let rmq = rmq::MaskRmq::new(&b.lcp);

    // Inclusive start pos in a of range of subtree, and right-lcp.
    let mut stack = vec![(0, 0)];
    for (i, &a_lcp_right) in a.lcp.iter().enumerate() {
        let mut start = i;
        while stack.last().unwrap().1 >= a_lcp_right {
            let (old_start, a_lcp) = stack.pop().unwrap();
            // Merge ranges.
            // TODO: use a faster algorithm. This is only nice if trees have depth O(log n).
            b_idx[old_start..=i].sort();

            // Update answer for all adjacent pairs from opposite ranges.
            // TODO: When ranges are very unbalanced, only check the insertion points.
            for (&(bl, al), &(br, ar)) in b_idx[old_start..=i].iter().tuple_windows() {
                if (al < start) ^ (ar < start) {
                    ans = max(ans, (a_lcp + rmq.query(bl..br), (al, ar)));
                }
            }

            start = old_start;
        }
        stack.push((start, a_lcp_right));
    }

    ans
}

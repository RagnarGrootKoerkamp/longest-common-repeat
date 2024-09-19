use std::cmp::min;

use coloured_trees::Tree;

use crate::{
    minimizers::{self},
    Ssa,
};

/// Find the longest common repeat length within a string.
pub fn lcr(t: &mut [u8], l: usize) -> usize {
    let k = min(l / 2, 64);
    let w = l - k - 1;
    // 1. Find minimizers
    let start = std::time::Instant::now();
    eprintln!("Find minimizers");

    let start = std::time::Instant::now();
    let mut minimizers = minimizers::minimizers_daniel(t, k, w);
    eprintln!("\tTime: {:?}", start.elapsed());
    eprintln!("Minimizers: {}", minimizers.len());
    eprintln!("1/density : {}", t.len() as f32 / minimizers.len() as f32);

    // eprintln!("Minimizers: {:?}", minimizers);

    // 2. SSA on minimizers.
    eprintln!("SSA1");
    let start = std::time::Instant::now();
    let ssa = Ssa::new(t, &minimizers);
    eprintln!("\tTime: {:?}", start.elapsed());
    // ssa.verify(t);
    // ssa.print(t);

    // 3. SSA on revere string minimizers.
    eprintln!("SSA2");
    let start = std::time::Instant::now();
    for i in &mut minimizers {
        *i = t.len() - *i;
    }
    t.reverse();
    let mut ssa_rev = Ssa::new(t, &minimizers);
    // ssa_rev.verify(t);
    // ssa_rev.print(t);
    t.reverse();
    for i in &mut ssa_rev.sa {
        *i = t.len() - *i;
    }
    eprintln!("\tTime: {:?}", start.elapsed());

    // Solve the common-tree problem.
    eprintln!("LCR");
    let start = std::time::Instant::now();
    let mut t1 = Tree {
        sa: ssa.sa,
        lcp: ssa.lcp,
    };
    t1.lcp.push(0);
    let mut t2 = Tree {
        sa: ssa_rev.sa,
        lcp: ssa_rev.lcp,
    };
    t2.lcp.push(0);

    // eprintln!("t1: {:?}", t1);
    // eprintln!("t2: {:?}", t2);

    // LCR.
    let (w, (n1, n2)) = coloured_trees::max_common_weight(&t1, &t2);
    eprintln!("\tTime: {:?}", start.elapsed());
    eprintln!("Weight: {w}, at positions {n1} and {n2}");
    w
}

#[cfg(test)]
mod test {
    use rand::{distributions::Alphanumeric, random, thread_rng, Rng};

    #[test]
    fn small() {
        let mut t = b"ABRACADABRAXYZPT".to_vec();
        let lcr = super::lcr(&mut t, 4);
        assert_eq!(lcr, 4);
    }
    #[test]
    fn large() {
        // Generate three random strings.
        let gen = |len: usize| {
            thread_rng()
                .sample_iter(&Alphanumeric)
                .take(len)
                .collect::<Vec<_>>()
        };
        let t1 = gen(random::<usize>() % 10);
        let t2 = gen(random::<usize>() % 10);
        let t3 = gen(random::<usize>() % 10);
        let secret = b"RagnarWasHere";
        let mut t = t1
            .iter()
            .chain(secret)
            .chain(&t2)
            .chain(secret)
            .chain(&t3)
            .cloned()
            .collect::<Vec<_>>();
        eprintln!("t: {}", t.len());

        // Lower bound on the length of the LCR we are looking for.
        let l = secret.len();
        let lcr = super::lcr(&mut t, l);
        assert!(
            lcr >= secret.len(),
            "LCR: {lcr} is not at least secret of length {}",
            secret.len()
        );
    }
}

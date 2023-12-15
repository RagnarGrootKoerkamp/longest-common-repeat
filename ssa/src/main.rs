use std::path::PathBuf;

use clap::Parser;
use rand::{seq::SliceRandom, thread_rng};
use ssa::Ssa;

#[derive(clap::Parser)]
struct Args {
    file: PathBuf,
    density: f32,
    l0: usize,
    #[clap(short, long)]
    exp_search: bool,
}

fn main() {
    let args = Args::parse();
    // Read std, drop lines starting with >, and concatenate the rest.
    eprintln!("Reading...");
    let t = std::fs::read_to_string(args.file)
        .unwrap()
        .lines()
        .skip(1)
        .collect::<String>();
    let t = t.as_bytes();
    eprintln!("Length: {}", t.len());
    // Read the density as the first command line argument.
    let density: f32 = args.density;
    // Take a subset of indices with the given density.
    let all_idxs = (0..t.len()).collect::<Vec<_>>();
    let idxs: Vec<_> = all_idxs
        .choose_multiple(&mut thread_rng(), (density * t.len() as f32) as usize)
        .cloned()
        .collect();
    eprintln!("Idxs: {}", idxs.len());
    eprintln!("Building... with exp_search {}", args.exp_search);
    let start = std::time::Instant::now();
    let _ssa = Ssa::new(t, idxs, Some(args.l0), args.exp_search);
    eprintln!("Time: {:?}", start.elapsed());
}

use std::{fmt::Write, path::PathBuf};

use clap::Parser;
use rand::{seq::SliceRandom, thread_rng};
use rdst::RadixSort;
use ssa::{lcr::lcr, Ssa};

#[derive(clap::Parser)]
struct Args {
    file: PathBuf,
    // density: f32,
    l0: usize,
    #[clap(short, long)]
    exp_search: bool,
}

fn main() {
    let args = Args::parse();
    // Read std, drop lines starting with >, and concatenate the rest.
    eprintln!("Reading...");
    let mut t = std::fs::read_to_string(&args.file)
        .unwrap()
        .lines()
        .filter(|s| !s.starts_with('>'))
        .map(|s| s.to_ascii_uppercase())
        .collect::<String>();
    let t = unsafe { t.as_bytes_mut() };
    eprintln!("Length: {}", t.len());

    lcr(t, args.l0);

    // Take a subset of indices with the given density.
    // let all_idxs = (0..t.len()).collect::<Vec<_>>();
    // let mut idxs: Vec<_> = all_idxs
    //     .choose_multiple(&mut thread_rng(), (args.density * t.len() as f32) as usize)
    //     .cloned()
    //     .collect();
    // idxs.sort_unstable();
    // let mut idxs_str = String::new();
    // for n in &idxs {
    //     writeln!(&mut idxs_str, "{}", n).unwrap();
    // }
    // std::fs::write(args.file.with_extension("idxs"), idxs_str).unwrap();
    // eprintln!("Idxs: {}", idxs.len());

    // eprintln!("Building... with exp_search {}", args.exp_search);
    // let start = std::time::Instant::now();
    // let ssa = Ssa::new_params(t, &idxs, Some(args.l0), args.exp_search);
    // eprintln!("Time: {:?}", start.elapsed());
    // eprintln!("Max LCS: {}", ssa.lcp.iter().max().unwrap());

    // eprintln!("Verifying...");
    // ssa.verify(t);

    // plot_lcp(&ssa);
}

fn plot_lcp(ssa: &Ssa) {
    let mut lcps = ssa.lcp.clone();
    lcps.radix_sort_unstable();
    let n = lcps.len();
    let max = lcps[n - 1];

    use plotters::prelude::*;
    let root = BitMapBackend::new("lcp.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .margin(5)
        .x_label_area_size(50)
        .y_label_area_size(50)
        .build_cartesian_2d(
            (1..max as u64).log_scale().zero_point(1 + max as u64),
            (0..n as u64).log_scale().zero_point(1 + n as u64),
        )
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            (0..n).map(|i| (lcps[i] as u64, i as u64)),
            RED,
        ))
        .unwrap();

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.9))
        .border_style(BLACK)
        .draw()
        .unwrap();

    root.present().unwrap();
}

//! This module contains benchmarks for polynomial multiplication
//! using both FFT-based and naive algorithms.
//!
//! Benchmarks are implemented using the `criterion` crate.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use fast_polynomial_arithmetic::{
    multiply_polynomials_fft, multiply_polynomials_naive, strategies,
};
use lambdaworks_math::fft::cpu::roots_of_unity::get_twiddles;
use lambdaworks_math::field::fields::fft_friendly::babybear_u32::Babybear31PrimeField;
use lambdaworks_math::field::traits::RootsConfig;
use proptest::prelude::*;
use proptest::strategy::ValueTree;
use proptest::test_runner::TestRunner;

// --- Combined Benchmarks ---

fn polynomial_multiplication_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Polynomial Multiplication Comparison");

    let mut runner = TestRunner::default();

    let degrees = vec![
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
        33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60,
        62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 110, 120,
        130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
        325, 350, 375, 400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
        1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
        2600, 2700, 2800, 2900, 3000,
    ];

    for &deg in degrees.iter() {
        // Benchmark FFT multiplication
        group.bench_with_input(BenchmarkId::new("FFT", deg), &deg, |b, &deg_val| {
            let strategy = (
                strategies::arb_polynomial(deg_val),
                strategies::arb_polynomial(deg_val),
            );
            b.iter_batched(
                || {
                    let (p1, p2) = strategy.new_tree(&mut runner).unwrap().current();
                    let min_domain_size = p1.degree() + p2.degree() + 1;
                    let n = strategies::next_power_of_2(min_domain_size);
                    let twiddles = get_twiddles::<Babybear31PrimeField>(
                        n.trailing_zeros() as u64,
                        RootsConfig::BitReverse,
                    )
                    .unwrap();
                    let inv_twiddles = get_twiddles::<Babybear31PrimeField>(
                        n.trailing_zeros() as u64,
                        RootsConfig::BitReverseInversed,
                    )
                    .unwrap();
                    (p1, p2, n, twiddles, inv_twiddles)
                },
                |(p1, p2, n, twiddles, inv_twiddles)| {
                    black_box(multiply_polynomials_fft(
                        &p1,
                        &p2,
                        n,
                        &twiddles,
                        &inv_twiddles,
                    ))
                },
                criterion::BatchSize::LargeInput,
            );
        });

        // Benchmark Naive multiplication
        group.bench_with_input(BenchmarkId::new("Naive", deg), &deg, |b, &deg_val| {
            let strategy = (
                strategies::arb_polynomial(deg_val),
                strategies::arb_polynomial(deg_val),
            );
            b.iter_batched(
                || {
                    let (p1, p2) = strategy.new_tree(&mut runner).unwrap().current();
                    (p1, p2)
                },
                |(p1, p2)| black_box(multiply_polynomials_naive(&p1, &p2)),
                criterion::BatchSize::LargeInput,
            );
        });
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = polynomial_multiplication_benchmark
}
criterion_main!(benches);

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
        10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450,
        500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 7500, 10000,
        12500, 15000, 17500, 20000,
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

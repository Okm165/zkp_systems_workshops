//! This crate provides functions for polynomial multiplication using
//! both the Fast Fourier Transform (FFT) algorithm and a naive O(N^2) approach.
//! It leverages the `lambdaworks_math` library for field arithmetic and FFT primitives.

use lambdaworks_math::fft::cpu::bit_reversing::in_place_bit_reverse_permute;
use lambdaworks_math::fft::cpu::fft::in_place_nr_2radix_fft;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::babybear_u32::Babybear31PrimeField;
use lambdaworks_math::polynomial::Polynomial;

// Type aliases for convenience, specifying the field to be Babybear31PrimeField.
type F = Babybear31PrimeField;
type FE = FieldElement<F>;

/// Multiplies two polynomials using the Fast Fourier Transform (FFT) algorithm.
///
/// This function performs polynomial multiplication in O(N log N) time, where N is
/// the size of the evaluation domain (a power of 2).
///
/// # Arguments
/// * `p1` - The first polynomial.
/// * `p2` - The second polynomial.
/// * `n` - The size of the FFT domain. Must be a power of 2 and `n >= degree(p1) + degree(p2) + 1`.
/// * `twiddles` - Precomputed bit-reversed roots of unity for the forward FFT.
/// * `inv_twiddles` - Precomputed bit-reversed inverse roots of unity for the Inverse FFT.
///
/// # Returns
/// A new `Polynomial` representing the product `p1 * p2`.
///
/// # Panics
/// This function does not explicitly panic, but relies on the correctness of `lambdaworks_math`
/// functions. Incorrect `n` or precomputed twiddles may lead to incorrect results.
pub fn multiply_polynomials_fft(
    p1: &Polynomial<FE>,
    p2: &Polynomial<FE>,
    n: usize,
    twiddles: &[FieldElement<F>],
    inv_twiddles: &[FieldElement<F>],
) -> Polynomial<FE> {
    // 1. Pad coefficients to match the FFT domain size `n`.
    // The FFT algorithm requires the input vectors to have a length equal to the domain size.
    let mut p1_coeffs = p1.coefficients.to_vec();
    p1_coeffs.resize(n, FE::zero()); // Pad with zeros.

    let mut p2_coeffs = p2.coefficients.to_vec();
    p2_coeffs.resize(n, FE::zero()); // Pad with zeros.

    // 2. Perform Fast Fourier Transform (FFT) on the padded coefficients.
    // The `in_place_nr_2radix_fft` function expects and produces bit-reversed evaluations
    // when used with bit-reversed twiddles.
    let mut p1_evals_bit_rev = p1_coeffs;
    in_place_nr_2radix_fft(&mut p1_evals_bit_rev, twiddles);

    let mut p2_evals_bit_rev = p2_coeffs;
    in_place_nr_2radix_fft(&mut p2_evals_bit_rev, twiddles);

    // 3. Apply bit-reversal permutation to get naturally ordered evaluations.
    // This step converts the bit-reversed output of the FFT into standard order.
    let mut p1_evals = p1_evals_bit_rev;
    in_place_bit_reverse_permute(&mut p1_evals);

    let mut p2_evals = p2_evals_bit_rev;
    in_place_bit_reverse_permute(&mut p2_evals);

    // 4. Perform pointwise multiplication of the evaluations.
    // This is the core step where the polynomial multiplication in the coefficient domain
    // is transformed into simple element-wise multiplication in the evaluation domain.
    let c_evals: Vec<FE> = p1_evals
        .iter()
        .zip(p2_evals.iter()) // Iterate over both evaluation vectors simultaneously.
        .map(|(y1, y2)| y1 * y2) // Multiply corresponding evaluations.
        .collect();

    // 5. Perform Inverse Fast Fourier Transform (IFFT) on the product evaluations.
    // This transforms the product evaluations back into product coefficients (still bit-reversed).
    let mut c_coeffs_bit_rev = c_evals;
    in_place_nr_2radix_fft(&mut c_coeffs_bit_rev, inv_twiddles);

    // 6. Apply bit-reversal permutation again to get naturally ordered, scaled coefficients.
    let mut c_coeffs_scaled = c_coeffs_bit_rev;
    in_place_bit_reverse_permute(&mut c_coeffs_scaled);

    // 7. Scale the coefficients by 1/N.
    // The IFFT process introduces a scaling factor of N (the domain size),
    // so we need to divide each coefficient by N to get the true coefficients.
    let n_inv = FE::from(n as u64)
        .inv()
        .expect("Inverse of N should exist in the field.");
    let c_coeffs: Vec<FE> = c_coeffs_scaled.iter().map(|c| c * n_inv).collect();

    // 8. Construct the resulting polynomial from the computed coefficients.
    Polynomial::new(&c_coeffs)
}

/// Multiplies two polynomials using a naive O(N^2) algorithm.
///
/// This function serves as a reference implementation for correctness verification
/// and to demonstrate the performance difference compared to the FFT-based approach.
///
/// # Arguments
/// * `p1` - The first polynomial.
/// * `p2` - The second polynomial.
///
/// # Returns
/// A new `Polynomial` representing the product `p1 * p2`.
pub fn multiply_polynomials_naive(p1: &Polynomial<FE>, p2: &Polynomial<FE>) -> Polynomial<FE> {
    let deg1 = p1.degree();
    let deg2 = p2.degree();

    // The degree of the product polynomial is deg1 + deg2.
    // The number of coefficients will be deg1 + deg2 + 1.
    let mut result_coeffs = vec![FE::zero(); deg1 + deg2 + 1];

    // Perform the standard polynomial multiplication by iterating through
    // each coefficient of p1 and multiplying it by each coefficient of p2.
    // The product of x^i and x^j contributes to the x^(i+j) term.
    for i in 0..=deg1 {
        for j in 0..=deg2 {
            result_coeffs[i + j] += p1.coefficients[i] * p2.coefficients[j];
        }
    }

    Polynomial::new(&result_coeffs)
}

pub mod strategies {
    use proptest::collection::vec;
    use proptest::prelude::{any, Strategy};

    use super::*;

    pub fn next_power_of_2(n: usize) -> usize {
        let mut k = 1;
        while k < n {
            k <<= 1;
        }
        k
    }

    /// Generates a polynomial with coefficients as `FE` elements,
    /// with a degree up to `max_degree`.
    pub fn arb_polynomial(max_degree: usize) -> impl Strategy<Value = Polynomial<FE>> {
        // Generate a vector of coefficients. The range `1..=max_degree`
        // ensures that the polynomial has at least one term (a constant).
        vec(any::<u64>().prop_map(FE::from), 1..=max_degree)
            .prop_map(|coeffs| Polynomial::new(&coeffs))
    }
}

#[cfg(test)]
mod tests {
    use lambdaworks_math::fft::cpu::roots_of_unity::get_twiddles;
    use lambdaworks_math::field::fields::fft_friendly::babybear_u32::Babybear31PrimeField;
    use lambdaworks_math::field::traits::RootsConfig;
    use proptest::prop_assert_eq;
    use proptest::test_runner::{Config, TestRunner};

    use crate::strategies::{self, next_power_of_2};
    use crate::{multiply_polynomials_fft, multiply_polynomials_naive};

    /// This test verifies that the FFT multiplication produces the same result
    /// as the naive multiplication for a range of randomly generated polynomials.
    #[test]
    fn proptest_fft_vs_naive_multiplication() {
        let mut runner = TestRunner::new(Config::default());

        let max_degree_for_proptest = 1000;

        let strategy = (
            strategies::arb_polynomial(max_degree_for_proptest),
            strategies::arb_polynomial(max_degree_for_proptest),
        );

        runner
            .run(&strategy, |(p1, p2)| {
                // Calculate expected result using the naive method.
                let expected_poly = multiply_polynomials_naive(&p1, &p2);

                // Determine the FFT domain size and precompute twiddles.
                let min_domain_size = p1.degree() + p2.degree() + 1;
                let n = next_power_of_2(min_domain_size);
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

                // Calculate actual result using the FFT method.
                let actual_poly = multiply_polynomials_fft(&p1, &p2, n, &twiddles, &inv_twiddles);

                // Assert that the coefficients are equal.
                prop_assert_eq!(
                    actual_poly.coefficients,
                    expected_poly.coefficients,
                    "FFT and Naive multiplication results differ!"
                );
                Ok(())
            })
            .unwrap();
    }
}

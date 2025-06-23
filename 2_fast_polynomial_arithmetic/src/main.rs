//! This `main` module demonstrates the usage of the FFT-based polynomial multiplication
//! with a concrete example using Babybear31PrimeField.
//! It sets up two polynomials, calculates their product using the `multiply_polynomials_fft`
//! function, and then verifies the result against a known expected polynomial.

// Import the FFT multiplication function from your library.
// Make sure this path is correct based on your crate structure.
use fast_polynomial_arithmetic::multiply_polynomials_fft;
use lambdaworks_math::fft::cpu::roots_of_unity::get_twiddles;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::babybear_u32::Babybear31PrimeField;
use lambdaworks_math::field::traits::RootsConfig;
use lambdaworks_math::polynomial::Polynomial;

// Type aliases for convenience, specifying the field.
type F = Babybear31PrimeField;
type FE = FieldElement<F>;

fn main() {
    // --- SETUP: DEFINE POLYNOMIALS AND DOMAIN ---

    // Define the first polynomial P1(x) = 3x^2 + 2x + 1
    let p1 = Polynomial::new(&[FE::from(1), FE::from(2), FE::from(3)]);
    // Define the second polynomial P2(x) = 5x + 2
    let p2 = Polynomial::new(&[FE::from(2), FE::from(5)]);

    // Calculate the minimum required domain size for FFT.
    // The degree of P1 is 2, and P2 is 1. The product C(x) will have degree 2 + 1 = 3.
    // We need a domain size N such that N >= degree(C) + 1. So, N >= 3 + 1 = 4.
    // The smallest power of 2 that satisfies this is N = 4.
    let n: usize = 4;

    println!("Multiplying polynomials using low-level FFT API:");
    println!("P1(x) = {}", p1.print_as_sage_poly(None));
    println!("P2(x) = {}", p2.print_as_sage_poly(None));
    println!("--------------------------------------\n");

    // Generate bit-reversed twiddle factors (roots of unity) for the forward FFT.
    // The `n.trailing_zeros()` gives log2(N), which is the 'k' for 2^k = N.
    let twiddles = get_twiddles::<F>(n.trailing_zeros() as u64, RootsConfig::BitReverse).unwrap();

    // Generate bit-reversed inverse twiddle factors for the Inverse FFT (IFFT).
    let inv_twiddles =
        get_twiddles::<F>(n.trailing_zeros() as u64, RootsConfig::BitReverseInversed).unwrap();

    // Perform the polynomial multiplication using the FFT algorithm.
    let c_poly = multiply_polynomials_fft(&p1, &p2, n, &twiddles, &inv_twiddles);

    // --- VERIFICATION ---
    println!("--- Verification ---");
    // Expected result: (3x^2 + 2x + 1) * (5x + 2) = 15x^3 + 16x^2 + 9x + 2
    // Coefficients are stored in ascending order of power: [constant, x^1, x^2, x^3]
    let expected_poly = Polynomial::new(&[FE::from(2), FE::from(9), FE::from(16), FE::from(15)]);

    println!(
        "Expected result: {}",
        expected_poly.print_as_sage_poly(None)
    );
    println!("Actual result:   {}", c_poly.print_as_sage_poly(None));

    // Assert that the coefficients of the computed polynomial match the expected coefficients.
    assert_eq!(
        c_poly.coefficients, expected_poly.coefficients,
        "The computed polynomial does not match the expected one."
    );
    println!("\nSuccess! The low-level FFT process produced the correct polynomial.");
}

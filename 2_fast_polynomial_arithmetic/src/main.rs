use lambdaworks_math::{
    fft::cpu::{
        bit_reversing::in_place_bit_reverse_permute, fft::in_place_nr_2radix_fft,
        roots_of_unity::get_twiddles,
    },
    field::{
        element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
        traits::RootsConfig,
    },
    polynomial::Polynomial,
};

type F = Stark252PrimeField;
type FE = FieldElement<F>;

fn main() {
    // --- SETUP: DEFINE POLYNOMIALS AND DOMAIN ---
    let p1 = Polynomial::new(&[FE::from(1), FE::from(2), FE::from(3)]); // P1(x) = 3x^2 + 2x + 1
    let p2 = Polynomial::new(&[FE::from(2), FE::from(5)]); // P2(x) = 5x + 2

    // As determined before, we need a domain of size N=4.
    let n = 4;
    // The order of the roots of unity is log2(N).
    let order = n;

    println!("Multiplying polynomials using low-level FFT API:");
    println!("P1(x) = {:?}", p1.print_as_sage_poly(None));
    println!("P2(x) = {:?}", p2.print_as_sage_poly(None));
    println!("--------------------------------------\n");

    // --- STEP 1: EVALUATION (Forward FFT) ---
    // This part explicitly demonstrates the DFT as a change of basis.

    // 1a. Pad coefficients to match the domain size.
    let mut p1_coeffs = p1.coefficients.to_vec();
    p1_coeffs.resize(n, FE::zero());
    let mut p2_coeffs = p2.coefficients.to_vec();
    p2_coeffs.resize(n, FE::zero());
    println!(
        "Step 1a: Padded P1 coeffs: {:?}",
        p1_coeffs.iter().map(|f| f.to_hex_str()).collect::<Vec<_>>()
    );
    println!(
        "         Padded P2 coeffs: {:?}\n",
        p2_coeffs.iter().map(|f| f.to_hex_str()).collect::<Vec<_>>()
    );

    // 1b. Generate the twiddle factors (the roots of unity).
    // The `in_place_nr_2radix_fft` function expects twiddles in bit-reversed order.
    let twiddles = get_twiddles::<F>(order as u64, RootsConfig::BitReverse).unwrap();
    println!("Step 1b: Generated bit-reversed twiddle factors.\n");

    // 1c. Perform the FFT.
    // The `in_place_nr_2radix_fft` function performs a "Natural to Reverse" transform.
    // It takes naturally ordered coefficients and outputs bit-reversed evaluations.
    let mut p1_evals_bit_rev = p1_coeffs.clone();
    in_place_nr_2radix_fft(&mut p1_evals_bit_rev, &twiddles);

    let mut p2_evals_bit_rev = p2_coeffs.clone();
    in_place_nr_2radix_fft(&mut p2_evals_bit_rev, &twiddles);
    println!("Step 1c: Performed core FFT transform. Output is in bit-reversed order.");
    println!(
        "         P1 bit-reversed evals: {:?}",
        p1_evals_bit_rev
            .iter()
            .map(|f| f.to_hex_str())
            .collect::<Vec<_>>()
    );
    println!(
        "         P2 bit-reversed evals: {:?}\n",
        p2_evals_bit_rev
            .iter()
            .map(|f| f.to_hex_str())
            .collect::<Vec<_>>()
    );

    // 1d. Apply the bit-reversal permutation to the output.
    // This reorders the evaluations into the natural order, making them ready for pointwise multiplication.
    // This is the explicit application of the permutation discussed in your text.
    let mut p1_evals = p1_evals_bit_rev;
    in_place_bit_reverse_permute(&mut p1_evals);

    let mut p2_evals = p2_evals_bit_rev;
    in_place_bit_reverse_permute(&mut p2_evals);
    println!("Step 1d: Applied bit-reversal permutation to get naturally ordered evaluations.");
    println!(
        "         P1 evaluations (natural order): {:?}",
        p1_evals.iter().map(|f| f.to_hex_str()).collect::<Vec<_>>()
    );
    println!(
        "         P2 evaluations (natural order): {:?}\n",
        p2_evals.iter().map(|f| f.to_hex_str()).collect::<Vec<_>>()
    );

    // --- STEP 2: POINTWISE MULTIPLICATION ---
    println!("Step 2: Performing pointwise multiplication of the evaluations.");
    let c_evals: Vec<FE> = p1_evals
        .iter()
        .zip(p2_evals.iter())
        .map(|(y1, y2)| y1 * y2)
        .collect();
    println!(
        "         C(x) evaluations (natural order): {:?}\n",
        c_evals.iter().map(|f| f.to_hex_str()).collect::<Vec<_>>()
    );

    // --- STEP 3: INTERPOLATION (Inverse FFT) ---
    // This mirrors the forward pass but uses inverse twiddles and includes a final scaling step.

    // 3a. Generate inverse twiddle factors.
    // As derived in your text, the IFFT uses ω⁻¹ instead of ω.
    let inv_twiddles = get_twiddles::<F>(order as u64, RootsConfig::BitReverseInversed).unwrap();
    println!("Step 3a: Generated bit-reversed inverse twiddle factors.\n");

    // 3b. Perform the Inverse FFT.
    // The structure is identical to the forward FFT.
    let mut c_coeffs_bit_rev = c_evals.clone();
    in_place_nr_2radix_fft(&mut c_coeffs_bit_rev, &inv_twiddles);
    println!("Step 3b: Performed core IFFT. Output coefficients are in bit-reversed order.");
    println!(
        "         C(x) bit-reversed coeffs: {:?}\n",
        c_coeffs_bit_rev
            .iter()
            .map(|f| f.to_hex_str())
            .collect::<Vec<_>>()
    );

    // 3c. Apply bit-reversal permutation to get coefficients in natural order.
    let mut c_coeffs_scaled = c_coeffs_bit_rev;
    in_place_bit_reverse_permute(&mut c_coeffs_scaled);
    println!("Step 3c: Applied bit-reversal permutation. Coefficients are now scaled.");
    println!(
        "         C(x) scaled coeffs: {:?}\n",
        c_coeffs_scaled
            .iter()
            .map(|f| f.to_hex_str())
            .collect::<Vec<_>>()
    );

    // 3d. Scale the result by 1/N.
    // This is the final step prescribed by the Inverse DFT theorem.
    let n_inv = FE::from(n as u64).inv().unwrap();
    let c_coeffs: Vec<FE> = c_coeffs_scaled.iter().map(|c| c * n_inv).collect();
    let c_poly = Polynomial::new(&c_coeffs);
    println!(
        "Step 3d: Scaled coefficients by 1/N = {}.",
        n_inv.to_hex_str()
    );
    println!(
        "         Final C(x) coeffs: {:?}\n",
        c_coeffs.iter().map(|f| f.to_hex_str()).collect::<Vec<_>>()
    );

    // --- VERIFICATION ---
    println!("--- Verification ---");
    let expected_poly = Polynomial::new(&[FE::from(2), FE::from(9), FE::from(16), FE::from(15)]);
    println!(
        "Expected result: {:?}",
        expected_poly.print_as_sage_poly(None)
    );
    println!("Actual result:   {:?}", c_poly.print_as_sage_poly(None));

    assert_eq!(c_poly.coefficients, expected_poly.coefficients);
    println!("Success! The low-level FFT process produced the correct polynomial.");
}

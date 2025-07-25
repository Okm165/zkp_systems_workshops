use lambdaworks_math::fft::cpu::roots_of_unity::get_powers_of_primitive_root_coset;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::babybear_u32::Babybear31PrimeField;
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::polynomial::Polynomial;

/// The prime field for our computations (Babybear).
type F = Babybear31PrimeField;
/// A field element in the Babybear field.
type FE = FieldElement<F>;

// ================================================================================================
// PART 1: COMPUTATION & EXECUTION TRACE
// ================================================================================================
// A STARK proof starts with a computation. We represent this computation as an "execution trace".
// For this example, our computation is the Fibonacci sequence.

/// Generates a Fibonacci sequence trace of a given length.
fn generate_fibonacci_trace(trace_length: usize) -> Vec<FE> {
    let mut trace = vec![FE::zero(); trace_length];
    // Set the initial values, which act as our boundary conditions.
    trace[0] = FE::one();
    trace[1] = FE::one();

    for i in 2..trace_length {
        trace[i] = trace[i - 1] + trace[i - 2];
    }
    trace
}

// ================================================================================================
// PART 2: ARITHMETIZATION
// ================================================================================================
// Arithmetization is the process of converting the execution trace into a set of polynomial
// constraints. If the constraints hold, the computation was performed correctly.

/// Holds the polynomials and domains related to the arithmetized trace.
struct Arithmetization {
    trace_length: usize,
    domain: Vec<FE>,
    domain_generator: FE,
    trace_poly: Polynomial<FE>,
    // The evaluations of the constraint polynomials over the LDE domain.
    // We store the evaluations directly to avoid interpolating and then re-evaluating,
    // which is more efficient.
    boundary_constraint_poly_lde: Vec<FE>,
    transition_constraint_poly_lde: Vec<FE>,
    // The domain used for low-degree extension (LDE).
    lde_domain: Vec<FE>,
}

impl Arithmetization {
    /// Performs the arithmetization of the execution trace.
    fn new(trace: &[FE], blowup_factor: usize) -> Self {
        println!("\n-- STEP 2: ARITHMETIZATION --------------------------------------");
        println!("The Prover transforms the execution trace into polynomial constraints.");

        let trace_length = trace.len();
        assert!(
            trace_length.is_power_of_two(),
            "Trace length must be a power of two for FFT-based interpolation."
        );

        // 1. Define the evaluation domain D_TRACE = {g^0, g^1, ..., g^(n-1)}.
        // This domain corresponds to the steps of our computation.
        let root_order = trace_length.trailing_zeros();
        let domain_generator = F::get_primitive_root_of_unity(root_order as u64).unwrap();
        let domain: Vec<FE> = (0..trace_length).map(|i| domain_generator.pow(i)).collect();

        // 2. Interpolate the trace over D_TRACE to get the trace polynomial t(x).
        // This creates a single polynomial whose evaluations at the domain points match the trace.
        // i.e., t(g^i) = trace[i] for i in [0, n-1].
        let trace_poly = Polynomial::interpolate_fft::<F>(trace).unwrap();
        println!(
            "  [2.1] Interpolated trace of {} elements into trace polynomial t(x) of degree {}.",
            trace_length,
            trace_poly.degree()
        );

        // 3. Define the LDE (Low-Degree Extension) Domain.
        // We evaluate our polynomials on a much larger domain to prevent a dishonest prover
        // from creating a fake polynomial that matches the constraints only on the small domain.
        let lde_domain_size = trace_length * blowup_factor;
        let lde_root_order = lde_domain_size.trailing_zeros();
        let lde_domain = get_powers_of_primitive_root_coset(
            lde_root_order as u64,
            lde_domain_size,
            &FE::from(3), /* A coset offset prevents zeroifiers evaluations equal to 0 (this
                           * would result in division by 0). */
        )
        .unwrap();

        // 4. Evaluate the trace polynomial on the LDE domain.
        // These evaluations, t_lde = {t(x) | x ∈ LDE_domain}, are what the Prover commits to.
        let trace_poly_lde = trace_poly.evaluate_slice(&lde_domain);

        // 5. Boundary Constraints: Ensure the computation starts and ends correctly.
        // Constraint: t(x) must be 1 at the first two steps (g^0 and g^1).
        // Polynomial form: B(x) = (t(x) - I(x)) / Z_B(x), where:
        // - I(x) is a polynomial that evaluates to 1 at g^0 and g^1.
        // - Z_B(x) = (x - g^0)(x - g^1) is a zerofier polynomial.
        // B(x) will be a polynomial (i.e., division is clean) iff the constraints hold.
        println!("  [2.2] Evaluating boundary constraints on the LDE domain...");
        let boundary_constraint_poly_lde = {
            let boundary_interpolant =
                Polynomial::interpolate(&[domain[0], domain[1]], &[FE::one(), FE::one()]).unwrap();
            let boundary_zerofier_poly = Polynomial::new(&[-domain[0], FE::one()])
                * Polynomial::new(&[-domain[1], FE::one()]);

            let numerator_lde = trace_poly_lde
                .iter()
                .zip(&lde_domain)
                .map(|(t_eval, x)| t_eval - boundary_interpolant.evaluate(x))
                .collect::<Vec<_>>();
            let denominator_lde = boundary_zerofier_poly.evaluate_slice(&lde_domain);

            let mut denominator_inv_lde = denominator_lde;
            FE::inplace_batch_inverse(&mut denominator_inv_lde).unwrap();

            numerator_lde
                .iter()
                .zip(denominator_inv_lde.iter())
                .map(|(num, den_inv)| num * den_inv)
                .collect::<Vec<_>>()
        };

        // 6. Transition Constraints: Ensure each step correctly follows from the previous ones.
        // Constraint: For Fibonacci, t(g^2 * x) = t(g * x) + t(x).
        // This must hold for all steps except the last two (where the next state is undefined).
        // Polynomial form: T(x) = (t(g^2 * x) - t(g * x) - t(x)) / Z_T(x), where:
        // - The numerator is the Fibonacci relation.
        // - Z_T(x) = (x^n - 1) / ((x - g^{n-2})(x - g^{n-1})) is the zerofier.
        // T(x) will be a polynomial iff the transition is valid for every step.
        println!("  [2.3] Evaluating transition constraints on the LDE domain...");
        let transition_constraint_poly_lde = {
            let trace_lde_g = trace_poly.evaluate_slice(
                &lde_domain
                    .iter()
                    .map(|x| x * domain_generator)
                    .collect::<Vec<_>>(),
            );
            let trace_lde_g2 = trace_poly.evaluate_slice(
                &lde_domain
                    .iter()
                    .map(|x| x * domain_generator.square())
                    .collect::<Vec<_>>(),
            );
            let numerator_lde = trace_lde_g2
                .iter()
                .zip(trace_lde_g.iter())
                .zip(trace_poly_lde.iter())
                .map(|((t_g2, t_g), t)| t_g2 - t_g - t)
                .collect::<Vec<_>>();

            // The zerofier Z_T(x) vanishes on all points of the trace domain except the
            // last two, where the transition constraint isn't supposed to hold.
            let transition_exemptions_poly =
                (Polynomial::new(&[-domain[trace_length - 2], FE::one()]))
                    * (Polynomial::new(&[-domain[trace_length - 1], FE::one()]));

            let mut exemptions_inv_lde = transition_exemptions_poly.evaluate_slice(&lde_domain);
            FE::inplace_batch_inverse(&mut exemptions_inv_lde).unwrap();

            // Z_T(x) = (x^n - 1) * Z_exemptions(x)^-1
            let denominator_lde = lde_domain
                .iter()
                .zip(exemptions_inv_lde.iter())
                .map(|(x, inv_exemption)| (x.pow(trace_length) - FE::one()) * inv_exemption)
                .collect::<Vec<_>>();

            let mut denominator_inv_lde = denominator_lde;
            FE::inplace_batch_inverse(&mut denominator_inv_lde).unwrap();

            numerator_lde
                .iter()
                .zip(denominator_inv_lde.iter())
                .map(|(num, den_inv)| num * den_inv)
                .collect::<Vec<_>>()
        };

        Self {
            trace_length,
            domain,
            domain_generator,
            trace_poly,
            boundary_constraint_poly_lde,
            transition_constraint_poly_lde,
            lde_domain,
        }
    }
}

// ================================================================================================
// PART 3: COMPOSITION & OUT-OF-DOMAIN SAMPLING
// ================================================================================================
// To reduce the number of polynomials the Verifier needs to check, the Prover combines them
// into a single "composition polynomial" H(x).

/// Holds the composition polynomial and its related data.
struct Composition {
    // LDE of the composition polynomial H(x).
    composition_poly_lde: Vec<FE>,
    // Coefficient form of H(x), needed for out-of-domain evaluation.
    composition_poly: Polynomial<FE>,
}

impl Composition {
    /// Combines constraint polynomials into a single composition polynomial H(x).
    /// H(x) = α₁ * B(x) + α₂ * T(x)
    /// The Verifier provides random challenges (α₁, α₂) to ensure the Prover can't cheat.
    fn new(arithmetization: &Arithmetization, alpha1: &FE, alpha2: &FE) -> Self {
        println!("\n-- STEP 3: POLYNOMIAL COMPOSITION -----------------------------");
        let composition_poly_lde = arithmetization
            .boundary_constraint_poly_lde
            .iter()
            .zip(&arithmetization.transition_constraint_poly_lde)
            .map(|(b_eval, t_eval)| b_eval * alpha1 + t_eval * alpha2)
            .collect::<Vec<_>>();

        // Interpolate to get the coefficient form. This is needed for the OOD check.
        let composition_poly =
            Polynomial::interpolate_offset_fft::<F>(&composition_poly_lde, &FE::from(3)).unwrap();

        println!(
            "  [3.1] Combined constraints into composition polynomial H(x) of degree {}.",
            composition_poly.degree()
        );
        println!("        The Prover commits to H(x) (e.g., via a Merkle tree of its LDE).");

        Self {
            composition_poly_lde,
            composition_poly,
        }
    }

    /// Simulates the out-of-domain check (the "DEEP" part of STARKs begins here).
    /// The Verifier asks the Prover to evaluate polynomials at a random point 'z' that is
    /// *not* in the LDE domain. This forces the Prover to have committed to actual low-degree
    /// polynomials, not just arbitrary values.
    fn perform_ood_check(
        &self,
        arithmetization: &Arithmetization,
        alpha1: &FE,
        alpha2: &FE,
        z: &FE,
    ) {
        println!("\n-- STEP 4: OUT-OF-DOMAIN SAMPLING (OOD) -----------------------");
        let g = &arithmetization.domain_generator;

        // Prover evaluates the trace polynomial at z and its required shifts (z*g, z*g^2),
        // and the composition polynomial H(z). These evaluations are sent to the verifier.
        let t_z = arithmetization.trace_poly.evaluate(z);
        let t_zg = arithmetization.trace_poly.evaluate(&(z * g));
        let t_zg2 = arithmetization.trace_poly.evaluate(&(z * g.square()));
        let h_z = self.composition_poly.evaluate(z);
        println!("  --> Prover to Verifier: Send evaluations at random point z.");
        println!(
            "      t(z)={}, t(z*g)={}, t(z*g^2)={}, H(z)={}",
            t_z.representative(),
            t_zg.representative(),
            t_zg2.representative(),
            h_z.representative()
        );

        // Verifier uses these evaluations to reconstruct H(z) on its own.
        // It computes the boundary and transition constraints at 'z' using the claimed t(z) values.
        println!("  <-- Verifier: Reconstructs H(z) to check consistency.");
        let boundary_interpolant = Polynomial::interpolate(
            &[arithmetization.domain[0], arithmetization.domain[1]],
            &[FE::one(), FE::one()],
        )
        .unwrap();
        let boundary_zerofier_z = (z - arithmetization.domain[0]) * (z - arithmetization.domain[1]);
        let boundary_eval_z =
            (t_z - boundary_interpolant.evaluate(z)) * boundary_zerofier_z.inv().unwrap();

        let transition_zerofier_z = {
            let numerator = z.pow(arithmetization.trace_length) - FE::one();
            let exemptions_at_z = (z - arithmetization.domain[arithmetization.trace_length - 2])
                * (z - arithmetization.domain[arithmetization.trace_length - 1]);
            numerator * exemptions_at_z.inv().unwrap()
        };
        let transition_eval_z = (t_zg2 - t_zg - t_z) * transition_zerofier_z.inv().unwrap();

        let h_z_reconstructed = boundary_eval_z * alpha1 + transition_eval_z * alpha2;

        println!(
            "      Reconstructed H(z): {}",
            h_z_reconstructed.representative()
        );
        assert_eq!(h_z, h_z_reconstructed, "Out-of-domain check failed!");
        println!("  [4.1] SUCCESS: Verifier's reconstructed H(z) matches Prover's H(z).");
    }
}

// ================================================================================================
// PART 4: DEEP COMPOSITION POLYNOMIAL
// ================================================================================================
// All the evaluation claims from the OOD check are bundled together into a single polynomial,
// the "DEEP" composition polynomial. Proving this single polynomial has a low degree is
// equivalent to proving all the original claims simultaneously.

/// Holds the DEEP composition polynomial.
struct DeepComposition {
    deep_poly_lde: Vec<FE>,
}

impl DeepComposition {
    /// Constructs the DEEP composition polynomial's evaluations over the LDE domain.
    /// D(x) = β₀ * (H(x) - H(z))/(x - z) +
    ///        β₁ * (t(x) - t(z))/(x - z) +
    ///        β₂ * (t(x) - t(zg))/(x - zg) +
    ///        β₃ * (t(x) - t(zg^2))/(x - zg^2)
    /// Each term will be a polynomial if and only if the numerator is zero when the
    /// denominator is zero (i.e., the evaluation claims are correct).
    fn new(
        arithmetization: &Arithmetization,
        composition: &Composition,
        z: &FE,
        betas: &[FE; 4],
    ) -> Self {
        println!("\n-- STEP 5: DEEP COMPOSITION -----------------------------------");
        println!("The Prover creates the DEEP polynomial D(x) to bundle all OOD claims.");
        let g = &arithmetization.domain_generator;

        // Pre-evaluate polynomials at OOD points (Prover already has these).
        let h_z = composition.composition_poly.evaluate(z);
        let t_z = arithmetization.trace_poly.evaluate(z);
        let t_zg = arithmetization.trace_poly.evaluate(&(z * g));
        let t_zg2 = arithmetization.trace_poly.evaluate(&(z * g.square()));

        // Get evaluations on LDE domain.
        let h_lde = &composition.composition_poly_lde;
        let t_lde = arithmetization
            .trace_poly
            .evaluate_slice(&arithmetization.lde_domain);

        // Compute point-wise evaluations for each term of the DEEP polynomial.
        let h_term_lde = h_lde
            .iter()
            .zip(&arithmetization.lde_domain)
            .map(|(h_xi, xi)| (h_xi - h_z) * (xi - z).inv().unwrap())
            .collect::<Vec<_>>();

        let t_z_term_lde = t_lde
            .iter()
            .zip(&arithmetization.lde_domain)
            .map(|(t_xi, xi)| (t_xi - t_z) * (xi - z).inv().unwrap())
            .collect::<Vec<_>>();

        let t_zg_term_lde = t_lde
            .iter()
            .zip(&arithmetization.lde_domain)
            .map(|(t_xi, xi)| (t_xi - t_zg) * (xi - (z * g)).inv().unwrap())
            .collect::<Vec<_>>();

        let t_zg2_term_lde = t_lde
            .iter()
            .zip(&arithmetization.lde_domain)
            .map(|(t_xi, xi)| (t_xi - t_zg2) * (xi - (z * g.square())).inv().unwrap())
            .collect::<Vec<_>>();

        // Combine the terms with random weights (betas) from the Verifier.
        let deep_poly_lde = h_term_lde
            .iter()
            .zip(&t_z_term_lde)
            .zip(&t_zg_term_lde)
            .zip(&t_zg2_term_lde)
            .map(|(((h_term, t_z_term), t_zg_term), t_zg2_term)| {
                h_term * betas[0]
                    + t_z_term * betas[1]
                    + t_zg_term * betas[2]
                    + t_zg2_term * betas[3]
            })
            .collect::<Vec<_>>();

        // The Prover would now run FRI on the `deep_poly_lde` to prove it has a low degree.
        // We will simulate the final check of that protocol.
        let deep_poly_coeffs =
            Polynomial::interpolate_offset_fft::<F>(&deep_poly_lde, &FE::from(3)).unwrap();
        println!(
            "  [5.1] Constructed DEEP polynomial D(x) of degree {}.",
            deep_poly_coeffs.degree()
        );
        println!(
            "        The Prover commits to D(x) and generates a FRI proof of its low-degreeness."
        );

        Self { deep_poly_lde }
    }

    /// Simulates the final spot-check after the FRI protocol.
    /// The FRI protocol gives the Verifier a random point `x₀` from the LDE domain and the
    /// claimed evaluations of the committed polynomials at that point. The Verifier checks if
    /// these values are consistent.
    fn perform_final_spot_check(
        &self,
        arithmetization: &Arithmetization,
        composition: &Composition,
        z: &FE,
        betas: &[FE; 4],
        x0_index: usize, // Index of a point in the LDE domain from a FRI query.
    ) {
        println!("\n-- STEP 6: FINAL CONSISTENCY CHECK --------------------------");
        let x0 = &arithmetization.lde_domain[x0_index];
        println!(
            "The Verifier picks a random point x₀={} from the LDE domain (via FRI).",
            x0.representative()
        );

        // Prover provides evaluations D(x₀), H(x₀), and t(x₀), authenticated by Merkle paths.
        let deep_x0 = self.deep_poly_lde[x0_index];
        let h_x0 = composition.composition_poly_lde[x0_index];
        let t_x0 = arithmetization.trace_poly.evaluate(x0);
        println!("  --> Prover to Verifier: Openings at x₀.");
        println!(
            "      D(x₀)={}, H(x₀)={}, t(x₀)={}",
            deep_x0.representative(),
            h_x0.representative(),
            t_x0.representative()
        );

        // Verifier reconstructs D(x₀) using the provided H(x₀), t(x₀) and the OOD
        // values it received earlier.
        println!("  <-- Verifier: Reconstructs D(x₀) to check final consistency.");
        let g = &arithmetization.domain_generator;

        // These OOD values are already known and trusted by the verifier from Step 4.
        let h_z = composition.composition_poly.evaluate(z);
        let t_z = arithmetization.trace_poly.evaluate(z);
        let t_zg = arithmetization.trace_poly.evaluate(&(z * g));
        let t_zg2 = arithmetization.trace_poly.evaluate(&(z * g.square()));

        let h_term_recon = (h_x0 - h_z) * (x0 - z).inv().unwrap();
        let t_z_term_recon = (t_x0 - t_z) * (x0 - z).inv().unwrap();
        let t_zg_term_recon = (t_x0 - t_zg) * (x0 - (z * g)).inv().unwrap();
        let t_zg2_term_recon = (t_x0 - t_zg2) * (x0 - (z * g.square())).inv().unwrap();

        let deep_x0_reconstructed = h_term_recon * betas[0]
            + t_z_term_recon * betas[1]
            + t_zg_term_recon * betas[2]
            + t_zg2_term_recon * betas[3];

        println!(
            "      Reconstructed D(x₀): {}",
            deep_x0_reconstructed.representative()
        );
        assert_eq!(deep_x0, deep_x0_reconstructed, "Final spot check failed!");
        println!("  [6.1] SUCCESS: All polynomial commitments are consistent.");
    }
}

/// This demo walks through the main algebraic steps of a STARK proving system,
/// from the initial computation to the final consistency checks. It omits the
/// cryptographic commitment scheme (Merkle Trees) and the FRI protocol itself,
/// focusing instead on the design of the polynomial constraints.
fn main() {
    println!("--- STARK Polynomial IOP Demo: Fibonacci Sequence ---");

    // ============================================================================
    // 1. SETUP & TRACE GENERATION
    // ============================================================================
    println!("\n-- STEP 1: EXECUTION TRACE ------------------------------------");
    let trace_length = 8;
    let blowup_factor = 8;
    let fib_trace = generate_fibonacci_trace(trace_length);
    println!("The Prover executes the computation and records the trace.");
    println!(
        "  Fibonacci trace (len {}): {:?}",
        trace_length,
        fib_trace
            .iter()
            .map(|e| e.representative())
            .collect::<Vec<_>>()
    );

    // ============================================================================
    // 2. PROVER: ARITHMETIZATION
    // ============================================================================
    let arithmetization = Arithmetization::new(&fib_trace, blowup_factor);

    // ============================================================================
    // 3. PROVER->VERIFIER INTERACTION
    // ============================================================================
    // Verifier sends random challenges to the Prover to ensure security.

    // <-- Verifier sends challenges α₁, α₂ for the composition polynomial.
    let alpha1 = FE::from(5);
    let alpha2 = FE::from(7);
    println!(
        "\n<-- Verifier to Prover: Send challenges α₁={}, α₂={}",
        alpha1.representative(),
        alpha2.representative()
    );
    let composition = Composition::new(&arithmetization, &alpha1, &alpha2);

    // <-- Verifier sends a random out-of-domain point 'z'.
    let z = FE::from(10);
    println!(
        "\n<-- Verifier to Prover: Send out-of-domain point z={}",
        z.representative()
    );
    composition.perform_ood_check(&arithmetization, &alpha1, &alpha2, &z);

    // <-- Verifier sends challenges β's for the DEEP polynomial.
    let betas = [FE::from(11), FE::from(13), FE::from(15), FE::from(17)];
    println!("\n<-- Verifier to Prover: Send challenges β's for DEEP polynomial");
    let deep_composition = DeepComposition::new(&arithmetization, &composition, &z, &betas);

    // ============================================================================
    // 4. FINAL VERIFICATION (Simulating a FRI query result)
    // ============================================================================
    // In a real system, the FRI protocol would conclude by querying a few points.
    // We simulate one such query and the final check.
    let x0_index = 5; // A random index into the LDE domain.
    deep_composition.perform_final_spot_check(&arithmetization, &composition, &z, &betas, x0_index);

    println!("\n\n--- Proof Verified Successfully ---");
}

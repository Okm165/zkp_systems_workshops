use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::babybear::Babybear31PrimeField;
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::polynomial::Polynomial;

/// The prime field for our computations (Babybear).
type F = Babybear31PrimeField;
/// A field element in the Babybear field.
type FE = FieldElement<F>;

// ================================================================================================
// PART 1: FROM COMPUTATION TO VERIFIABLE CLAIMS
// ================================================================================================
// We represent the Fibonacci sequence computation as an execution trace.

/// Generates a Fibonacci sequence trace of a given length.
fn generate_fibonacci_trace(trace_length: usize) -> Vec<FE> {
    let mut trace = vec![FE::zero(); trace_length];
    trace[0] = FE::one();
    trace[1] = FE::one();

    for i in 2..trace_length {
        trace[i] = &trace[i - 1] + &trace[i - 2];
    }
    trace
}

// ================================================================================================
// PART 2: ARITHMETIZATION AND CONSTRAINT DESIGN
// ================================================================================================
// We convert the trace into polynomials and define constraints algebraically.

struct Arithmetization {
    trace_length: usize,
    domain: Vec<FE>,
    trace_poly: Polynomial<FE>,
    // Store LDEs for use in composition.
    boundary_constraint_poly_lde: Vec<FE>,
    transition_constraint_poly_lde: Vec<FE>,
    // LDE domain for evaluations.
    lde_domain: Vec<FE>,
}

impl Arithmetization {
    /// Performs the arithmetization of the execution trace.
    fn new(trace: &[FE], blowup_factor: usize) -> Self {
        println!("- Step 2.1: Arithmetization -");
        let trace_length = trace.len();
        assert!(trace_length.is_power_of_two());

        // 1. Generate the evaluation domain D_S = {g^0, g^1, ..., g^{n-1}}
        let root_order = trace_length.trailing_zeros();
        let generator = F::get_primitive_root_of_unity(root_order as u64).unwrap();
        let domain: Vec<FE> = (0..trace_length).map(|i| generator.pow(i)).collect();

        // 2. Interpolate the trace to get the trace polynomial t(x)
        let trace_poly = Polynomial::interpolate_fft::<F>(trace).unwrap();
        println!(
            "  - Interpolated trace of {} elements into t(x) of degree {}.",
            trace_length,
            trace_poly.degree()
        );

        // 3. Define LDE (Low-Degree Extension) Domain
        let lde_domain_size = trace_length * blowup_factor;
        let lde_root_order = lde_domain_size.trailing_zeros();
        let lde_generator = F::get_primitive_root_of_unity(lde_root_order as u64).unwrap();
        let lde_domain: Vec<FE> = (0..lde_domain_size)
            .map(|i| lde_generator.pow(i) * FE::from(3))
            .collect();

        // 4. Evaluate Trace Polynomial on LDE Domain
        let trace_poly_lde = trace_poly.evaluate_slice(&lde_domain);

        // 5. Design and Evaluate Boundary Constraints
        let boundary_constraint_poly_lde = {
            let boundary_interpolant = Polynomial::interpolate(
                &[domain[0].to_owned(), domain[1].to_owned()],
                &[FE::one(), FE::one()],
            )
            .unwrap();
            let boundary_zerofier = Polynomial::new(&[-domain[0].to_owned(), FE::one()])
                * Polynomial::new(&[-domain[1].to_owned(), FE::one()]);

            // Evaluate numerator and denominator on LDE domain
            let numerator_lde = trace_poly_lde
                .iter()
                .zip(&lde_domain)
                .map(|(t_eval, x)| t_eval - boundary_interpolant.evaluate(x))
                .collect::<Vec<_>>();
            let denominator_lde = boundary_zerofier.evaluate_slice(&lde_domain);

            // Point-wise division
            let mut denominator_inv_lde = denominator_lde;
            FE::inplace_batch_inverse(&mut denominator_inv_lde).unwrap();

            let evals = numerator_lde
                .iter()
                .zip(denominator_inv_lde.iter())
                .map(|(num, den_inv)| num * den_inv)
                .collect::<Vec<_>>();

            evals
        };

        // 6. Design and Evaluate Transition Constraints
        let transition_constraint_poly_lde = {
            // Numerator: t(x * g^2) - t(x * g) - t(x)
            let trace_lde_g = trace_poly.evaluate_slice(
                &lde_domain
                    .iter()
                    .map(|x| x * &generator)
                    .collect::<Vec<_>>(),
            );
            let trace_lde_g2 = trace_poly.evaluate_slice(
                &lde_domain
                    .iter()
                    .map(|x| x * &generator.square())
                    .collect::<Vec<_>>(),
            );
            let numerator_lde = trace_lde_g2
                .iter()
                .zip(trace_lde_g.iter())
                .zip(trace_poly_lde.iter())
                .map(|((t_g2, t_g), t)| t_g2 - t_g - t)
                .collect::<Vec<_>>();

            // Denominator (Zerofier): Z_T(x) = (x^n - 1) / ((x - g^{n-2})(x-g^{n-1}))
            let zerofier_numerator = Polynomial::new_monomial(FE::one(), trace_length) - FE::one();
            let exemptions = (Polynomial::new(&[-domain[trace_length - 2].to_owned(), FE::one()]))
                * (Polynomial::new(&[-domain[trace_length - 1].to_owned(), FE::one()]));

            let mut exemptions_inv_lde = exemptions.evaluate_slice(&lde_domain);
            FE::inplace_batch_inverse(&mut exemptions_inv_lde).unwrap();

            let denominator_lde = zerofier_numerator
                .evaluate_slice(&lde_domain)
                .iter()
                .zip(exemptions_inv_lde.iter())
                .map(|(n, d_inv)| n * d_inv)
                .collect::<Vec<_>>();

            // Point-wise division
            let mut denominator_inv_lde = denominator_lde;
            FE::inplace_batch_inverse(&mut denominator_inv_lde).unwrap();

            let evals = numerator_lde
                .iter()
                .zip(denominator_inv_lde.iter())
                .map(|(num, den_inv)| num * den_inv)
                .collect::<Vec<_>>();

            evals
        };

        Self {
            trace_length,
            domain,
            trace_poly,
            boundary_constraint_poly_lde,
            transition_constraint_poly_lde,
            lde_domain,
        }
    }
}

// ================================================================================================
// PART 3: THE COMPOSITION POLYNOMIAL AND OUT-OF-DOMAIN CHECK
// ================================================================================================
// We combine all constraints and perform a crucial check to prevent spoofing.

struct Composition {
    // Store LDE for DEEP poly construction.
    composition_poly_lde: Vec<FE>,
    // Store coefficient form for OOD evaluation.
    composition_poly: Polynomial<FE>,
}

impl Composition {
    /// Combines constraint polynomials into a single composition polynomial.
    fn new(arithmetization: &Arithmetization, beta1: &FE, beta2: &FE) -> Self {
        println!("\n- Step 3.1: Composition -");

        // Construct the LDE of H(x) directly from the LDEs of B(x) and C(x).
        let composition_poly_lde = arithmetization
            .boundary_constraint_poly_lde
            .iter()
            .zip(arithmetization.transition_constraint_poly_lde.iter())
            .map(|(b, c)| b * beta1 + c * beta2)
            .collect::<Vec<_>>();

        // Interpolate to get the coefficient form.
        let composition_poly = Polynomial::interpolate_fft::<F>(&composition_poly_lde).unwrap();

        println!(
            "  - Constructed composition polynomial H(x) of degree {}.",
            composition_poly.degree()
        );

        Self {
            composition_poly_lde,
            composition_poly,
        }
    }

    /// Simulates the out-of-domain check.
    fn perform_ood_check(&self, arithmetization: &Arithmetization, beta1: &FE, beta2: &FE, z: &FE) {
        println!("\n- Step 3.2: Out-of-Domain Check -");
        let generator =
            F::get_primitive_root_of_unity(arithmetization.trace_length.trailing_zeros() as u64)
                .unwrap();

        // Prover provides evaluations at z and its shifts.
        let t_z = arithmetization.trace_poly.evaluate(z);
        let t_zg = arithmetization.trace_poly.evaluate(&(z * &generator));
        let t_zg2 = arithmetization
            .trace_poly
            .evaluate(&(z * &generator.square()));
        let h_z = self.composition_poly.evaluate(z);
        println!("  - Prover sends evaluations at random point z: t(z), t(zg), t(zg^2), H(z)");

        // Verifier reconstructs H(z) from trace evaluations.
        let boundary_interpolant = Polynomial::interpolate(
            &[
                arithmetization.domain[0].to_owned(),
                arithmetization.domain[1].to_owned(),
            ],
            &[FE::one(), FE::one()],
        )
        .unwrap();
        let boundary_zerofier_z =
            (z - arithmetization.domain[0].to_owned()) * (z - arithmetization.domain[1].to_owned());
        let boundary_eval =
            (&t_z - boundary_interpolant.evaluate(z)) * boundary_zerofier_z.inv().unwrap();

        let transition_zerofier_z = {
            let numerator = z.pow(arithmetization.trace_length) - FE::one();
            let denominator = (z - arithmetization.domain[arithmetization.trace_length - 2]
                .to_owned())
                * (z - arithmetization.domain[arithmetization.trace_length - 1].to_owned());
            numerator * denominator.inv().unwrap()
        };
        let transition_eval = (&t_zg2 - &t_zg - &t_z) * transition_zerofier_z.inv().unwrap();

        let h_z_reconstructed = &boundary_eval * beta1 + &transition_eval * beta2;

        println!("  - Verifier reconstructs H(z) using trace evaluations.");
        assert_eq!(h_z, h_z_reconstructed, "Out-of-domain check failed!");
        println!("  - SUCCESS: Verifier's reconstructed H(z) matches Prover's H(z).");
    }
}

// ================================================================================================
// PART 4: THE DEEP COMPOSITION POLYNOMIAL AND FINAL VERIFICATION
// ================================================================================================
// We batch all evaluation proofs into a single polynomial for a final FRI proof.

struct DeepComposition {
    deep_poly_lde: Vec<FE>,
}

impl DeepComposition {
    /// Constructs the DEEP composition polynomial's evaluations over the LDE domain.
    fn new(
        arithmetization: &Arithmetization,
        composition: &Composition,
        z: &FE,
        gammas: &[FE; 4],
    ) -> Self {
        println!("\n- Step 4.1: DEEP Composition Polynomial Construction -");
        let g =
            F::get_primitive_root_of_unity(arithmetization.trace_length.trailing_zeros() as u64)
                .unwrap();

        let lde_domain = &arithmetization.lde_domain;

        // Pre-evaluate polynomials at OOD points
        let h_z = composition.composition_poly.evaluate(z);
        let t_z = arithmetization.trace_poly.evaluate(z);
        let t_zg = arithmetization.trace_poly.evaluate(&(z * &g));
        let t_zg2 = arithmetization.trace_poly.evaluate(&(z * &g.square()));

        // Get evaluations on LDE domain
        let h_lde = &composition.composition_poly_lde;
        let t_lde = arithmetization.trace_poly.evaluate_slice(lde_domain);

        // Compute point-wise evaluations for each term of the DEEP polynomial
        let h_term_lde = lde_domain
            .iter()
            .zip(h_lde.iter())
            .map(|(x_i, h_x_i)| (h_x_i - &h_z) * (x_i - z).inv().unwrap() * &gammas[0])
            .collect::<Vec<_>>();

        let t_z_term_lde = lde_domain
            .iter()
            .zip(t_lde.iter())
            .map(|(x_i, t_x_i)| (t_x_i - &t_z) * (x_i - z).inv().unwrap() * &gammas[1])
            .collect::<Vec<_>>();

        let t_zg_term_lde = lde_domain
            .iter()
            .zip(t_lde.iter())
            .map(|(x_i, t_x_i)| (t_x_i - &t_zg) * (x_i - (z * &g)).inv().unwrap() * &gammas[2])
            .collect::<Vec<_>>();

        let t_zg2_term_lde = lde_domain
            .iter()
            .zip(t_lde.iter())
            .map(|(x_i, t_x_i)| {
                (t_x_i - &t_zg2) * (x_i - (z * &g.square())).inv().unwrap() * &gammas[3]
            })
            .collect::<Vec<_>>();

        // Sum the terms to get the final DEEP LDE
        let deep_poly_lde = h_term_lde
            .iter()
            .zip(t_z_term_lde.iter())
            .zip(t_zg_term_lde.iter())
            .zip(t_zg2_term_lde.iter())
            .map(|(((h, t_z), t_zg), t_zg2)| h + t_z + t_zg + t_zg2)
            .collect::<Vec<_>>();

        Self { deep_poly_lde }
    }

    /// Simulates the final spot-check.
    fn perform_final_spot_check(
        &self,
        arithmetization: &Arithmetization,
        composition: &Composition,
        z: &FE,
        gammas: &[FE; 4],
        x0_index: usize, // In-domain point index from FRI query
    ) {
        println!("\n- Step 4.2: Final Spot Check at in-domain point x0 -");

        let x0 = &arithmetization.lde_domain[x0_index];

        // Prover provides evaluations
        let deep_x0 = self.deep_poly_lde[x0_index].clone();
        let h_x0 = composition.composition_poly_lde[x0_index].clone();
        let t_x0 = arithmetization.trace_poly.evaluate(x0);
        println!("  - Prover opens Deep(x0), H(x0), and t(x0) with Merkle proofs.");

        // Verifier reconstructs Deep(x0)
        let g =
            F::get_primitive_root_of_unity(arithmetization.trace_length.trailing_zeros() as u64)
                .unwrap();

        // These values are already known by the verifier from the OOD check
        let h_z = composition.composition_poly.evaluate(z);
        let t_z = arithmetization.trace_poly.evaluate(z);
        let t_zg = arithmetization.trace_poly.evaluate(&(z * &g));
        let t_zg2 = arithmetization.trace_poly.evaluate(&(z * &g.square()));

        let h_term_recon = (&h_x0 - &h_z) * (x0 - z).inv().unwrap();
        let t_z_term_recon = (&t_x0 - &t_z) * (x0 - z).inv().unwrap();
        let t_zg_term_recon = (&t_x0 - &t_zg) * (x0 - (z * &g)).inv().unwrap();
        let t_zg2_term_recon = (&t_x0 - &t_zg2) * (x0 - (z * &g.square())).inv().unwrap();

        let deep_x0_reconstructed = &h_term_recon * &gammas[0]
            + &t_z_term_recon * &gammas[1]
            + &t_zg_term_recon * &gammas[2]
            + &t_zg2_term_recon * &gammas[3];

        println!("  - Verifier reconstructs Deep(x0) from trusted values.");
        assert_eq!(deep_x0, deep_x0_reconstructed, "Final spot check failed!");
        println!("  - SUCCESS: Verifier's reconstructed Deep(x0) matches Prover's Deep(x0).");
    }
}

fn main() {
    println!("--- Running AIR Constraint Design Demo ---");

    // ============================================================================
    // PART 1: SETUP
    // ============================================================================
    let trace_length = 8; // Use a small power of 2 for this example
    let blowup_factor = 8;
    let fib_trace = generate_fibonacci_trace(trace_length);
    println!("\n- Step 1: Execution Trace Generation -");
    println!(
        "  - Fibonacci trace (length {}): {:?}",
        trace_length,
        fib_trace
            .iter()
            .map(|e| e.representative())
            .collect::<Vec<_>>()
    );

    // ============================================================================
    // PART 2: ARITHMETIZATION
    // ============================================================================
    let arithmetization = Arithmetization::new(&fib_trace, blowup_factor);

    // ============================================================================
    // PART 3: COMPOSITION & OOD CHECK
    // ============================================================================
    // Verifier sends challenges β₁, β₂
    let beta1 = FE::from(5);
    let beta2 = FE::from(7);
    println!(
        "\n<-- Verifier sends challenges β₁={}, β₂={}",
        beta1.representative(),
        beta2.representative()
    );

    let composition = Composition::new(&arithmetization, &beta1, &beta2);

    // Verifier sends random out-of-domain point z
    let z = FE::from(10);
    println!(
        "\n<-- Verifier sends out-of-domain point z={}",
        z.representative()
    );

    composition.perform_ood_check(&arithmetization, &beta1, &beta2, &z);

    // ============================================================================
    // PART 4: DEEP COMPOSITION & FINAL CHECK
    // ============================================================================
    // Verifier sends challenges γ₁, γ₂, γ₃, γ₄
    let gammas = [FE::from(11), FE::from(13), FE::from(15), FE::from(17)];
    println!("\n<-- Verifier sends challenges γ's for DEEP polynomial");

    let deep_composition = DeepComposition::new(&arithmetization, &composition, &z, &gammas);

    // Verifier chooses a random in-domain point x₀ (from FRI queries)
    // We'll simulate this by picking an index into the LDE domain.
    let x0_index = 5;
    println!(
        "\n<-- Verifier chooses in-domain spot-check point x₀={}",
        arithmetization.lde_domain[x0_index].representative()
    );

    deep_composition.perform_final_spot_check(
        &arithmetization,
        &composition,
        &z,
        &gammas,
        x0_index,
    );

    println!("\n--- Demo Finished Successfully ---");
}

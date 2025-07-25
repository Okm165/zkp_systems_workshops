// ================================================================================================
// PART 2: ARITHMETIZATION
// ================================================================================================
// Arithmetization is the process of converting the execution trace into a set of polynomial
// constraints. If the constraints hold, the computation was performed correctly.

use lambdaworks_math::fft::cpu::roots_of_unity::get_powers_of_primitive_root_coset;
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::polynomial::Polynomial;

use crate::{F, FE};

/// Holds the polynomials and domains related to the arithmetized trace.
pub struct Arithmetization {
    pub trace_length: usize,
    pub domain: Vec<FE>,
    pub domain_generator: FE,
    pub trace_poly: Polynomial<FE>,
    // The evaluations of the constraint polynomials over the LDE domain.
    // We store the evaluations directly to avoid interpolating and then re-evaluating,
    // which is more efficient.
    pub boundary_constraint_poly_lde: Vec<FE>,
    pub transition_constraint_poly_lde: Vec<FE>,
    // The domain used for low-degree extension (LDE).
    pub lde_domain: Vec<FE>,
}

impl Arithmetization {
    /// Performs the arithmetization of the execution trace.
    pub fn new(trace: &[FE], blowup_factor: usize) -> Self {
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

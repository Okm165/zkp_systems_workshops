// ================================================================================================
// PART 3: COMPOSITION & OUT-OF-DOMAIN SAMPLING
// ================================================================================================
// To reduce the number of polynomials the Verifier needs to check, the Prover combines them
// into a single "composition polynomial" H(x).

use lambdaworks_math::polynomial::Polynomial;

use crate::arithmetization::Arithmetization;
use crate::{F, FE};

/// Holds the composition polynomial and its related data.
pub struct Composition {
    // LDE of the composition polynomial H(x).
    pub composition_poly_lde: Vec<FE>,
    // Coefficient form of H(x), needed for out-of-domain evaluation.
    pub composition_poly: Polynomial<FE>,
}

impl Composition {
    /// Combines constraint polynomials into a single composition polynomial H(x).
    /// H(x) = α₁ * B(x) + α₂ * T(x)
    /// The Verifier provides random challenges (α₁, α₂) to ensure the Prover can't cheat.
    pub fn new(arithmetization: &Arithmetization, alpha1: &FE, alpha2: &FE) -> Self {
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
    pub fn perform_ood_check(
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

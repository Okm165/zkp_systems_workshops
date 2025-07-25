// ================================================================================================
// PART 4: DEEP COMPOSITION POLYNOMIAL
// ================================================================================================
// All the evaluation claims from the OOD check are bundled together into a single polynomial,
// the "DEEP" composition polynomial. Proving this single polynomial has a low degree is
// equivalent to proving all the original claims simultaneously.

use lambdaworks_math::polynomial::Polynomial;

use crate::arithmetization::Arithmetization;
use crate::composition::Composition;
use crate::{F, FE};

/// Holds the DEEP composition polynomial.
pub struct DeepComposition {
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
    pub fn new(
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
    pub fn perform_final_spot_check(
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

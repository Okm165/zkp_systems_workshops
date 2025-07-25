use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::babybear_u32::Babybear31PrimeField;

use crate::arithmetization::Arithmetization;
use crate::composition::Composition;
use crate::deep_composition::DeepComposition;
use crate::trace::generate_fibonacci_trace;

pub mod arithmetization;
pub mod composition;
pub mod deep_composition;
pub mod trace;

/// The prime field for our computations (Babybear).
type F = Babybear31PrimeField;
/// A field element in the Babybear field.
type FE = FieldElement<F>;

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

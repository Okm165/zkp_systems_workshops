//! # Educational FRI Protocol Implementation
//!
//! This code provides a simplified, educational implementation of the FRI (Fast Reed-Solomon
//! Interactive Oracle Proof of Proximity) protocol in Rust. It is designed for teaching
//! purposes to demonstrate the core concepts of FRI, which is a foundational component in
//! many modern STARK (Scalable Transparent Argument of Knowledge) systems.
//!
//! The implementation uses the `lambdaworks` library for finite field arithmetic, polynomials,
//! and Merkle trees.
//!
//! ## Protocol Flow Overview
//!
//! 1. **COMMIT**: The Prover evaluates a polynomial `P(x)` over a large domain (a Low-Degree
//!    Extension or LDE). It then commits to these evaluations using a Merkle tree.
//!
//! 2. **FOLD**: The Prover and Verifier engage in a recursive process. In each round:
//!     - The Verifier sends a random challenge, `beta`.
//!     - The Prover uses `beta` to "fold" the current set of evaluations into a smaller set,
//!       representing a new polynomial of half the degree.
//!     - The Prover commits to the new evaluations and the process repeats.
//!
//! 3. **LAST LAYER**: This folding continues until the polynomial is reduced to a constant. The
//!    Prover sends this constant value to the Verifier.
//!
//! 4. **QUERY**: The Verifier asks the Prover to reveal the evaluations of the polynomial at
//!    specific random points from the initial domain, along with their Merkle authentication paths
//!    for all layers.
//!
//! 5. **VERIFY**: The Verifier checks two things:
//!     - **Merkle Paths**: That the revealed evaluations are consistent with the commitments.
//!     - **Folding Consistency**: That the folding process was performed correctly at each step for
//!       the queried points. This ensures the Prover didn't cheat during the folding phase.
use lambdaworks_crypto::merkle_tree::backends::types::Keccak256Backend;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::babybear::Babybear31PrimeField;
use lambdaworks_math::polynomial::Polynomial;

use crate::prover::Prover;
use crate::types::FriParameters;
use crate::verifier::Verifier;

pub mod error;
pub mod prover;
pub mod types;
pub mod verifier;

/// The prime field for our computations (Babybear).
type F = Babybear31PrimeField;
/// A field element in the Babybear field.
type FE = FieldElement<F>;
/// The backend for our Merkle Tree, using Keccak256 for hashing.
type FriBackend = Keccak256Backend<F>;
/// The name of the protocol, used for initializing the transcript.
const PROTOCOL_ID: &[u8] = b"Educational FRI";

fn main() {
    // 1. SETUP
    // The polynomial we want to prove knowledge of: P(x) = x^3 - 3x + 2
    let poly = Polynomial::new(&[FE::from(2), -FE::from(3), FE::from(0), FE::from(1)]);
    let claimed_degree = 3;
    // Parameters: degree 3, blowup factor 4 (domain size 16), 2 queries.
    let params = FriParameters::new(claimed_degree, 8, 2);

    // 2. PROVE
    // The Prover generates a proof that it knows a polynomial of degree <= 3.
    let mut prover = Prover::new(poly, params.clone());
    let proof = prover.prove().unwrap();

    // 3. VERIFY
    // The Verifier checks the proof.
    let mut verifier = Verifier::new(params);
    match verifier.verify(&proof) {
        Ok(_) => println!("\n✅ SUCCESS: Proof verified successfully!"),
        Err(e) => println!("\n❌ FAILURE: Proof verification failed: {}", e),
    }
}

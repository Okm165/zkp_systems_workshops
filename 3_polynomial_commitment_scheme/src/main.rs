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

// --- IMPORTS ---
use std::fmt;

use lambdaworks_crypto::fiat_shamir::default_transcript::DefaultTranscript;
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_crypto::merkle_tree::backends::types::Keccak256Backend;
use lambdaworks_crypto::merkle_tree::merkle::MerkleTree;
use lambdaworks_crypto::merkle_tree::proof::Proof;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::babybear::Babybear31PrimeField;
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::polynomial::Polynomial;
use lambdaworks_math::traits::AsBytes;

// --- TYPE ALIASES & CONSTANTS ---

/// The prime field for our computations (Babybear).
type F = Babybear31PrimeField;
/// A field element in the Babybear field.
type FE = FieldElement<F>;
/// The backend for our Merkle Tree, using Keccak256 for hashing.
type FriBackend = Keccak256Backend<F>;
/// The name of the protocol, used for initializing the transcript.
const PROTOCOL_ID: &[u8] = b"Educational FRI";

// --- CUSTOM ERROR TYPE ---

/// Defines the possible errors that can occur during the FRI protocol execution.
#[derive(Debug, PartialEq, Eq)]
pub enum FriError {
    /// Error building a Merkle tree for a specific layer.
    MerkleTreeConstructionError(String),
    /// A Merkle proof verification failed.
    InvalidMerkleProof,
    /// The folding process at a specific layer was inconsistent.
    InconsistentFolding {
        layer: usize,
        expected: String,
        got: String,
    },
}

impl fmt::Display for FriError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FriError::MerkleTreeConstructionError(msg) => {
                write!(f, "Merkle tree construction failed: {}", msg)
            }
            FriError::InvalidMerkleProof => write!(f, "Invalid Merkle proof"),
            FriError::InconsistentFolding {
                layer,
                expected,
                got,
            } => write!(
                f,
                "Inconsistent folding at layer {}: expected {}, got {}",
                layer, expected, got
            ),
        }
    }
}

// --- SHARED PROTOCOL PARAMETERS & STRUCTURES ---

/// Shared parameters for the FRI protocol, agreed upon by the Prover and Verifier.
#[derive(Debug, Clone)]
pub struct FriParameters {
    /// The initial evaluation domain (LDE).
    pub domain_0: Vec<FE>,
    /// How many queries the Verifier will make to check the proof.
    pub num_queries: usize,
    /// The size of the initial evaluation domain.
    pub domain_0_size: usize,
}

impl FriParameters {
    /// Creates a new set of FRI parameters.
    ///
    /// # Arguments
    /// * `claimed_degree`: The claimed degree of the initial polynomial.
    /// * `blowup_factor`: How much larger the evaluation domain is than the number of coefficients.
    ///   A larger factor provides more security.
    /// * `num_queries`: The number of queries to perform. More queries also increase security.
    pub fn new(claimed_degree: usize, blowup_factor: usize, num_queries: usize) -> Self {
        // The Low-Degree Extension (LDE) domain size.
        let domain_size = (claimed_degree + 1) * blowup_factor;
        // The domain is a multiplicative subgroup, so its size must be a power of 2.
        let root_order = domain_size.trailing_zeros();
        let generator = F::get_primitive_root_of_unity(root_order as u64).unwrap();
        let domain_0 = (0..domain_size).map(|i| generator.pow(i)).collect();

        Self {
            domain_0,
            num_queries,
            domain_0_size: domain_size,
        }
    }
}

/// Represents a single layer in the FRI protocol's commitment-folding process.
#[derive(Clone)]
struct FriLayer {
    /// The evaluations of the polynomial for this layer.
    evaluations: Vec<FE>,
    /// The Merkle tree committing to the evaluations.
    merkle_tree: MerkleTree<FriBackend>,
    /// The domain over which the evaluations were made.
    domain: Vec<FE>,
}

/// A decommitment for a single query, providing evaluations and Merkle paths for each layer.
#[derive(Debug, Clone)]
struct QueryDecommitment {
    /// The evaluation at the query index `q` for each layer.
    layer_evaluations: Vec<FE>,
    /// The Merkle authentication path for `layer_evaluations` at each layer.
    layer_auth_paths: Vec<Vec<[u8; 32]>>,
    /// The evaluation at the symmetric index `-q` for each layer.
    layer_evaluations_sym: Vec<FE>,
    /// The Merkle authentication path for `layer_evaluations_sym` at each layer.
    layer_auth_paths_sym: Vec<Vec<[u8; 32]>>,
}

/// The complete FRI proof sent from the Prover to the Verifier.
#[derive(Debug, Clone)]
pub struct FriProof {
    /// The Merkle root of each FRI layer.
    layer_commitments: Vec<[u8; 32]>,
    /// The value of the final, constant polynomial.
    last_layer_value: FE,
    /// The decommitments for each query.
    query_decommitments: Vec<QueryDecommitment>,
}

// --- CORE FOLDING LOGIC ---

/// Folds a layer of evaluations based on a challenge `beta`.
/// This is the heart of the FRI protocol's recursive step.
///
/// It takes a polynomial `f(x)` represented by its evaluations over a domain `D`,
/// and computes the evaluations of a new, smaller polynomial `f_next(x^2)` over `D^2`.
///
/// The formula is: `f_next(x^2) = (f(x) + f(-x))/2 + beta * (f(x) - f(-x))/(2x)`
/// where `(f(x) + f(-x))/2` is the even part of `f` and `(f(x) - f(-x))/(2x)` is the odd part.
fn fold_evaluations(evaluations: &[FE], domain: &[FE], beta: &FE) -> (Vec<FE>, Vec<FE>) {
    let next_domain_size = domain.len() / 2;
    let two_inv = FE::from(2).inv().unwrap();

    let next_evaluations = (0..next_domain_size)
        .map(|i| {
            // Get the evaluation at a point x and its symmetric counterpart -x
            let y = &evaluations[i];
            let y_symmetric = &evaluations[i + next_domain_size]; // Corresponds to -x

            // Get the domain value x and its inverse
            let x = &domain[i];
            let x_inv = x.inv().unwrap();

            // Calculate the even and odd components of the polynomial
            let f_even = (y + y_symmetric) * &two_inv;
            let f_odd = (y - y_symmetric) * &two_inv * &x_inv;

            // Combine them to get the evaluation of the next polynomial
            f_even + beta * f_odd
        })
        .collect();

    // The next domain consists of the squares of the first half of the current domain
    let next_domain = domain
        .iter()
        .take(next_domain_size)
        .map(|x| x.square())
        .collect();

    (next_evaluations, next_domain)
}

// --- PROVER ---

/// The Prover entity for the FRI protocol.
pub struct Prover {
    poly: Polynomial<FE>,
    params: FriParameters,
    transcript: DefaultTranscript<F>,
}

impl Prover {
    /// Creates a new Prover.
    pub fn new(poly: Polynomial<FE>, params: FriParameters) -> Self {
        Self {
            poly,
            params,
            transcript: DefaultTranscript::new(PROTOCOL_ID),
        }
    }

    /// Executes the entire proving process.
    pub fn prove(&mut self) -> Result<FriProof, FriError> {
        println!("--- Prover: Starting proof generation ---");

        // 1. Commit Phase: Evaluate the polynomial and commit to the evaluations.
        let initial_layer = self.commit_phase()?;
        // 2. Fold Phase: Recursively fold the polynomial until it's a constant.
        let (layers, last_value) = self.fold_phase(initial_layer)?;
        // 3. Query Phase: Generate decommitments for random queries.
        let query_decommitments = self.query_phase(&layers);

        println!("--- Prover: Proof generation complete ---\n");
        Ok(FriProof {
            layer_commitments: layers.iter().map(|l| l.merkle_tree.root).collect(),
            last_layer_value: last_value,
            query_decommitments,
        })
    }

    /// Phase 1: Commit to the initial polynomial evaluations on the LDE domain.
    fn commit_phase(&mut self) -> Result<FriLayer, FriError> {
        println!("[Prover] Phase 1: COMMIT");
        // Evaluate the polynomial on the large domain (LDE).
        let evaluations = self.poly.evaluate_slice(&self.params.domain_0);
        // Build a Merkle tree from the evaluations to commit to them.
        let merkle_tree = MerkleTree::<FriBackend>::build(&evaluations).ok_or_else(|| {
            FriError::MerkleTreeConstructionError("Failed to build initial Merkle tree".to_string())
        })?;

        // Add the Merkle root to the transcript to make it part of the public record.
        self.transcript.append_bytes(&merkle_tree.root);
        println!(
            "  > Layer 0 committed with root: 0x{}",
            hex::encode(merkle_tree.root)
        );

        Ok(FriLayer {
            evaluations,
            merkle_tree,
            domain: self.params.domain_0.clone(),
        })
    }

    /// Phase 2: Interactively fold the polynomial evaluations until a constant is reached.
    fn fold_phase(&mut self, initial_layer: FriLayer) -> Result<(Vec<FriLayer>, FE), FriError> {
        println!("[Prover] Phase 2: FOLD");
        let mut layers = vec![initial_layer];

        // Continue folding until the polynomial becomes a constant (evaluations list has 1 element)
        while layers.last().unwrap().evaluations.len() > 1 {
            let i = layers.len() - 1;
            // Get a random challenge `beta` from the transcript.
            let beta: FE = self.transcript.sample_field_element();
            println!(
                "  > Round {}: Sampled challenge beta = {}",
                i,
                beta.representative()
            );

            let previous_layer = layers.last().unwrap();
            // Fold the evaluations and domain for the next layer.
            let (next_evaluations, next_domain) =
                fold_evaluations(&previous_layer.evaluations, &previous_layer.domain, &beta);

            // Commit to the new evaluations.
            let next_merkle_tree =
                MerkleTree::<FriBackend>::build(&next_evaluations).ok_or_else(|| {
                    FriError::MerkleTreeConstructionError(format!(
                        "Failed to build Merkle tree for layer {}",
                        i + 1
                    ))
                })?;

            // Add the new Merkle root to the transcript.
            self.transcript.append_bytes(&next_merkle_tree.root);
            println!(
                "    - Layer {} committed with root: 0x{}",
                i + 1,
                hex::encode(next_merkle_tree.root)
            );

            layers.push(FriLayer {
                evaluations: next_evaluations,
                merkle_tree: next_merkle_tree,
                domain: next_domain,
            });
        }

        // The final layer contains a single evaluation, which is the constant value.
        let last_value = layers.last().unwrap().evaluations[0].clone();
        self.transcript.append_bytes(&last_value.as_bytes());
        println!(
            "  > Folding complete. Final value: {}",
            last_value.representative()
        );

        Ok((layers, last_value))
    }

    /// Phase 3: Generate decommitments for random queries issued by the verifier.
    fn query_phase(&mut self, layers: &[FriLayer]) -> Vec<QueryDecommitment> {
        println!("[Prover] Phase 3: QUERY");
        // Sample random indices from the transcript for the queries.
        let query_indices: Vec<usize> = (0..self.params.num_queries)
            .map(|_| self.sample_index(self.params.domain_0_size))
            .collect();

        println!(
            "  > Generating decommitments for queries at indices: {:?}",
            query_indices
        );

        query_indices
            .into_iter()
            .map(|mut query_idx| {
                let mut decommitment = QueryDecommitment {
                    layer_evaluations: Vec::new(),
                    layer_auth_paths: Vec::new(),
                    layer_evaluations_sym: Vec::new(),
                    layer_auth_paths_sym: Vec::new(),
                };

                // For each layer, provide the evaluation and its Merkle proof.
                for layer in layers {
                    let domain_size = layer.domain.len();
                    // The symmetric index corresponds to f(-x).
                    let sym_idx = (query_idx + domain_size / 2) % domain_size;

                    // Provide evaluation and auth path for f(x).
                    decommitment
                        .layer_evaluations
                        .push(layer.evaluations[query_idx].clone());
                    decommitment.layer_auth_paths.push(
                        layer
                            .merkle_tree
                            .get_proof_by_pos(query_idx)
                            .unwrap()
                            .merkle_path,
                    );

                    // Provide evaluation and auth path for f(-x).
                    decommitment
                        .layer_evaluations_sym
                        .push(layer.evaluations[sym_idx].clone());
                    decommitment.layer_auth_paths_sym.push(
                        layer
                            .merkle_tree
                            .get_proof_by_pos(sym_idx)
                            .unwrap()
                            .merkle_path,
                    );

                    // The index for the next layer is `query_idx mod (domain_size / 2)`.
                    query_idx %= (domain_size / 2).max(1);
                }
                decommitment
            })
            .collect()
    }

    /// Samples a random index from the transcript.
    fn sample_index(&mut self, max_value: usize) -> usize {
        // Use 8 bytes from the transcript for a u64, then get a value in range.
        let sample_bytes: [u8; 8] = self.transcript.sample()[..8].try_into().unwrap();
        (u64::from_be_bytes(sample_bytes) % max_value as u64) as usize
    }
}

// --- VERIFIER ---

/// The Verifier entity for the FRI protocol.
pub struct Verifier {
    params: FriParameters,
    transcript: DefaultTranscript<F>,
}

impl Verifier {
    /// Creates a new Verifier.
    pub fn new(params: FriParameters) -> Self {
        Self {
            params,
            transcript: DefaultTranscript::new(PROTOCOL_ID),
        }
    }

    /// Verifies the FRI proof.
    pub fn verify(&mut self, proof: &FriProof) -> Result<(), FriError> {
        println!("--- Verifier: Starting verification ---");

        // Reconstruct the challenges (`betas`) and query indices by replaying the transcript.
        let (betas, query_indices) = self.reconstruct_challenges(proof);

        let root_order = self.params.domain_0_size.trailing_zeros();
        let generator = F::get_primitive_root_of_unity(root_order as u64).unwrap();

        // Verify each query independently.
        for (query_num, &query_idx) in query_indices.iter().enumerate() {
            println!(
                "\n[Verifier] Verifying query #{} (for original index {})",
                query_num + 1,
                query_idx
            );
            self.verify_query(
                proof,
                query_idx,
                &betas,
                &generator,
                &proof.query_decommitments[query_num],
            )?;
        }

        Ok(())
    }

    /// Reconstructs all challenges by replaying the Prover's commitments from the proof.
    /// This ensures the Verifier uses the exact same random values as the Prover.
    fn reconstruct_challenges(&mut self, proof: &FriProof) -> (Vec<FE>, Vec<usize>) {
        // Feed the commitments into the transcript in the same order as the Prover.
        self.transcript.append_bytes(&proof.layer_commitments[0]);
        let betas: Vec<FE> = proof
            .layer_commitments
            .iter()
            .skip(1)
            .map(|commitment| {
                // Sample the field element *before* appending the next commitment.
                let beta = self.transcript.sample_field_element();
                self.transcript.append_bytes(commitment);
                beta
            })
            .collect();

        // Feed the last layer's value.
        self.transcript
            .append_bytes(&proof.last_layer_value.as_bytes());

        // Now, sample the query indices. They will be the same as the Prover's.
        let query_indices = (0..proof.query_decommitments.len())
            .map(|_| self.sample_index(self.params.domain_0_size))
            .collect();

        println!("[Verifier] Reconstructed challenges and query indices from proof commitments.");
        (betas, query_indices)
    }

    /// Verifies a single query decommitment.
    fn verify_query(
        &self,
        proof: &FriProof,
        query_idx: usize,
        betas: &[FE],
        generator: &FE,
        decommitment: &QueryDecommitment,
    ) -> Result<(), FriError> {
        // Step 1: Verify the Merkle proofs for each layer's evaluations.
        self.verify_merkle_paths(proof, query_idx, decommitment)?;

        // Step 2: Verify the folding consistency across all layers.
        self.verify_folding_consistency(proof, query_idx, decommitment, betas, generator)?;

        Ok(())
    }

    /// Verifies that all evaluations in a decommitment are valid against the layer commitments.
    fn verify_merkle_paths(
        &self,
        proof: &FriProof,
        query_idx: usize,
        decommitment: &QueryDecommitment,
    ) -> Result<(), FriError> {
        let mut current_idx = query_idx;

        for i in 0..proof.layer_commitments.len() {
            let domain_size = self.params.domain_0_size >> i;
            let sym_idx = (current_idx + domain_size / 2) % domain_size;
            let commitment = &proof.layer_commitments[i];

            // Verify proof for f(x)
            let proof_path = Proof {
                merkle_path: decommitment.layer_auth_paths[i].clone(),
            };
            if !proof_path.verify::<FriBackend>(
                commitment,
                current_idx,
                &decommitment.layer_evaluations[i],
            ) {
                return Err(FriError::InvalidMerkleProof);
            }

            // Verify proof for f(-x)
            let proof_path_sym = Proof {
                merkle_path: decommitment.layer_auth_paths_sym[i].clone(),
            };
            if !proof_path_sym.verify::<FriBackend>(
                commitment,
                sym_idx,
                &decommitment.layer_evaluations_sym[i],
            ) {
                return Err(FriError::InvalidMerkleProof);
            }

            println!(
                "  > Layer {}: Merkle proofs valid for indices {} and {}",
                i, current_idx, sym_idx
            );
            current_idx %= (domain_size / 2).max(1);
        }
        Ok(())
    }

    /// Checks that the folding from layer `i` to `i+1` was done correctly.
    fn verify_folding_consistency(
        &self,
        proof: &FriProof,
        query_idx: usize,
        decommitment: &QueryDecommitment,
        betas: &[FE],
        generator: &FE,
    ) -> Result<(), FriError> {
        // Start with the claimed evaluation from the *next* layer and work backwards.
        // `claimed_child_evaluation` is the value at layer `i+1` that we are checking.
        let mut claimed_child_evaluation = proof.last_layer_value.clone();

        // Iterate backwards from the second-to-last layer down to the first.
        for i in (0..proof.layer_commitments.len() - 1).rev() {
            // Get the evaluations for f(x) and f(-x) at the current layer `i`.
            let y = &decommitment.layer_evaluations[i];
            let y_sym = &decommitment.layer_evaluations_sym[i];

            // Recompute `x` for the specific query index at this layer's domain size.
            let domain_size = self.params.domain_0_size >> i;
            let current_query_idx_in_layer = query_idx % domain_size;
            let g_i = generator.pow(1_u64 << i); // Generator for the i-th domain
            let x = g_i.pow(current_query_idx_in_layer);
            let x_inv = x.inv().unwrap();

            // Re-compute what the folded value should be using the folding formula.
            let two_inv = FE::from(2).inv().unwrap();
            let f_even = (y + y_sym) * &two_inv;
            let f_odd = (y - y_sym) * &two_inv * &x_inv;
            let expected_child_evaluation = &f_even + &betas[i] * &f_odd;

            // Check if our calculation matches the claimed evaluation from the next layer.
            if claimed_child_evaluation != expected_child_evaluation {
                return Err(FriError::InconsistentFolding {
                    layer: i,
                    expected: expected_child_evaluation.representative().to_hex(),
                    got: claimed_child_evaluation.representative().to_hex(),
                });
            }

            println!("  > Layer {}->{}: Folding is consistent.", i, i + 1);

            // For the next iteration, the "child" becomes the current evaluation.
            claimed_child_evaluation = y.clone();
        }

        Ok(())
    }

    /// Samples a random index from the transcript.
    fn sample_index(&mut self, max_value: usize) -> usize {
        // Use 8 bytes from the transcript for a u64, then get a value in range.
        let sample_bytes: [u8; 8] = self.transcript.sample()[..8].try_into().unwrap();
        (u64::from_be_bytes(sample_bytes) % max_value as u64) as usize
    }
}

// --- MAIN EXECUTION ---

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

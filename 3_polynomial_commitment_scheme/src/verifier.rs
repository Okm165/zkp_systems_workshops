use lambdaworks_crypto::fiat_shamir::default_transcript::DefaultTranscript;
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_crypto::merkle_tree::proof::Proof;
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::traits::AsBytes;

use crate::error::FriError;
use crate::types::{FriParameters, FriProof, QueryDecommitment};
use crate::{FriBackend, F, FE, PROTOCOL_ID};

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

        let root_order = self.params.domain.len().trailing_zeros();
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
            .map(|_| self.sample_index(self.params.domain.len()))
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
            let domain_size = self.params.domain.len() >> i;
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
            let domain_size = self.params.domain.len() >> i;
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

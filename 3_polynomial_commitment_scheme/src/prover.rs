use lambdaworks_crypto::fiat_shamir::default_transcript::DefaultTranscript;
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_crypto::merkle_tree::merkle::MerkleTree;
use lambdaworks_math::polynomial::Polynomial;
use lambdaworks_math::traits::AsBytes;

use crate::error::FriError;
use crate::types::{FriLayer, FriParameters, FriProof, QueryDecommitment};
use crate::{FriBackend, F, FE, PROTOCOL_ID};

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
                Self::fold_evaluations(&previous_layer.evaluations, &previous_layer.domain, &beta);

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

    /// Samples a random index from the transcript.
    fn sample_index(&mut self, max_value: usize) -> usize {
        // Use 8 bytes from the transcript for a u64, then get a value in range.
        let sample_bytes: [u8; 8] = self.transcript.sample()[..8].try_into().unwrap();
        (u64::from_be_bytes(sample_bytes) % max_value as u64) as usize
    }
}

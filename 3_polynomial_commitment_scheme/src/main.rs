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

// Type aliases for convenience, specifying the field to be Babybear31PrimeField.
type F = Babybear31PrimeField;
type FE = FieldElement<F>;
type FriBackend = Keccak256Backend<F>;

/// Shared parameters for the FRI protocol, agreed upon by the Prover and Verifier.
#[derive(Debug, Clone)]
struct FriParameters {
    domain_0: Vec<FE>,
    num_queries: usize,
}

impl FriParameters {
    /// Creates a new set of FRI parameters.
    pub fn new(claimed_degree: usize, blowup_factor: usize, num_queries: usize) -> Self {
        let domain_size = (claimed_degree + 1) * blowup_factor;
        let root_order = domain_size.trailing_zeros();
        let generator = F::get_primitive_root_of_unity(root_order as u64).unwrap();
        let domain_0 = (0..domain_size).map(|i| generator.pow(i)).collect();

        Self {
            domain_0,
            num_queries,
        }
    }

    pub fn domain_0_size(&self) -> usize {
        self.domain_0.len()
    }
}

/// Represents a single layer in the FRI protocol.
#[derive(Clone)]
struct FriLayer {
    evaluations: Vec<FE>,
    merkle_tree: MerkleTree<FriBackend>,
    domain: Vec<FE>,
}

/// A decommitment for a single query, providing evaluations and Merkle paths for each layer.
#[derive(Debug, Clone)]
struct QueryDecommitment {
    layer_evaluations: Vec<FE>,
    layer_auth_paths: Vec<Vec<[u8; 32]>>,
    layer_evaluations_sym: Vec<FE>,
    layer_auth_paths_sym: Vec<Vec<[u8; 32]>>,
}

/// The complete FRI proof sent from the Prover to the Verifier.
#[derive(Debug, Clone)]
struct FriProof {
    layer_commitments: Vec<[u8; 32]>,
    last_layer_value: FE,
    query_decommitments: Vec<QueryDecommitment>,
}

/// Folds a layer of evaluations based on a challenge `beta`.
/// This is the core recursive step of the FRI protocol.
fn fold_evaluations(evaluations: &[FE], domain: &[FE], beta: &FE) -> (Vec<FE>, Vec<FE>) {
    let next_domain_size = domain.len() / 2;
    let two_inv = FE::from(2).inv().unwrap();

    let next_evaluations = (0..next_domain_size)
        .map(|i| {
            let x = domain[i].clone();
            let x_inv = x.inv().unwrap();
            let y = &evaluations[i];
            let y_symmetric = &evaluations[i + next_domain_size]; // Corresponds to -x

            // f_next(x^2) = (f(x) + f(-x))/2 + beta * (f(x) - f(-x))/(2x)
            let f_even = (y + y_symmetric) * &two_inv;
            let f_odd = (y - y_symmetric) * &two_inv * &x_inv;

            f_even + beta * f_odd
        })
        .collect();

    let next_domain = domain
        .iter()
        .take(next_domain_size)
        .map(|x| x.square())
        .collect();

    (next_evaluations, next_domain)
}

/// The Prover entity for the FRI protocol.
struct Prover {
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
            transcript: DefaultTranscript::new(b"FRI example transcript"),
        }
    }

    /// Executes the entire proving process.
    pub fn prove(&mut self) -> Result<FriProof, String> {
        println!("--- Prover: Starting proof generation ---");

        let initial_layer = self.commit_phase()?;
        let (layers, last_value) = self.fold_phase(initial_layer)?;
        let query_decommitments = self.query_phase(&layers)?;

        Ok(FriProof {
            layer_commitments: layers.iter().map(|l| l.merkle_tree.root).collect(),
            last_layer_value: last_value,
            query_decommitments,
        })
    }

    /// Phase 1: Commit to the initial polynomial evaluations on the LDE domain.
    fn commit_phase(&mut self) -> Result<FriLayer, String> {
        println!("  Phase 1: COMMIT");
        let evaluations = self.poly.evaluate_slice(&self.params.domain_0);
        let merkle_tree = MerkleTree::<FriBackend>::build(&evaluations)
            .ok_or("Failed to build initial Merkle tree")?;

        self.transcript.append_bytes(&merkle_tree.root);
        println!(
            "    Layer 0: Committed with root {:x}",
            merkle_tree
                .root
                .iter()
                .fold(0, |acc, &x| (acc << 8) | x as u32)
        );

        Ok(FriLayer {
            evaluations,
            merkle_tree,
            domain: self.params.domain_0.clone(),
        })
    }

    /// Phase 2: Interactively fold the polynomial evaluations until a constant is reached.
    fn fold_phase(&mut self, initial_layer: FriLayer) -> Result<(Vec<FriLayer>, FE), String> {
        println!("  Phase 2: FOLD");
        let mut layers = vec![initial_layer];
        let log_domain_size = self.params.domain_0_size().trailing_zeros();

        for i in 0..log_domain_size {
            let beta: FE = self.transcript.sample_field_element();
            println!(
                "    Round {}: Sampled challenge beta = {}",
                i,
                beta.representative()
            );

            let previous_layer = layers.last().unwrap();
            let (next_evaluations, next_domain) =
                fold_evaluations(&previous_layer.evaluations, &previous_layer.domain, &beta);

            if next_domain.is_empty() {
                break;
            }

            let next_merkle_tree = MerkleTree::<FriBackend>::build(&next_evaluations)
                .ok_or_else(|| format!("Failed to build Merkle tree for layer {}", i + 1))?;

            self.transcript.append_bytes(&next_merkle_tree.root);
            println!(
                "      Layer {}: Committed with root {:x}",
                i + 1,
                next_merkle_tree
                    .root
                    .iter()
                    .fold(0, |acc, &x| (acc << 8) | x as u32)
            );

            layers.push(FriLayer {
                evaluations: next_evaluations,
                merkle_tree: next_merkle_tree,
                domain: next_domain,
            });
        }

        let last_value = layers.last().unwrap().evaluations[0].clone();
        self.transcript.append_bytes(&last_value.as_bytes());
        println!(
            "  Folding complete. Final value: {}",
            last_value.representative()
        );

        Ok((layers, last_value))
    }

    /// Phase 3: Generate decommitments for random queries issued by the verifier.
    fn query_phase(&mut self, layers: &[FriLayer]) -> Result<Vec<QueryDecommitment>, String> {
        println!("  Phase 3: QUERY");
        let query_indices: Vec<usize> = (0..self.params.num_queries)
            .map(|_| {
                (u64::from_be_bytes(self.transcript.sample()[..8].try_into().unwrap())
                    % self.params.domain_0_size() as u64) as usize
            })
            .collect();
        println!(
            "    Generating decommitments for queries at indices: {:?}",
            query_indices
        );

        Ok(query_indices
            .into_iter()
            .map(|mut query_idx| {
                let mut decommitment = QueryDecommitment {
                    layer_evaluations: Vec::new(),
                    layer_auth_paths: Vec::new(),
                    layer_evaluations_sym: Vec::new(),
                    layer_auth_paths_sym: Vec::new(),
                };

                for layer in layers {
                    let domain_size = layer.domain.len();
                    let sym_idx = (query_idx + domain_size / 2) % domain_size;

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

                    query_idx %= (domain_size / 2).max(1);
                }
                decommitment
            })
            .collect())
    }
}

/// The Verifier entity for the FRI protocol.
struct Verifier {
    params: FriParameters,
    transcript: DefaultTranscript<F>,
}

#[derive(Debug, PartialEq, Eq)]
enum FriValidationError {
    InvalidMerkleProof,
    InconsistentFolding { layer: usize, expected: String, got: String },
}

impl Verifier {
    /// Creates a new Verifier.
    pub fn new(params: FriParameters) -> Self {
        Self {
            params,
            transcript: DefaultTranscript::new(b"FRI example transcript"),
        }
    }

    /// Verifies the FRI proof.
    pub fn verify(&mut self, proof: &FriProof) -> Result<(), FriValidationError> {
        println!("\n--- Verifier: Starting verification ---");

        let (betas, query_indices) = self.reconstruct_challenges(proof);

        let root_order = self.params.domain_0_size().trailing_zeros();
        let generator = F::get_primitive_root_of_unity(root_order as u64).unwrap();
        for (query_num, &query_idx) in query_indices.iter().enumerate() {
            println!(
                "  Verifying query #{} for original index {}",
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
    fn reconstruct_challenges(&mut self, proof: &FriProof) -> (Vec<FE>, Vec<usize>) {
        self.transcript.append_bytes(&proof.layer_commitments[0]);
        let betas: Vec<FE> = proof
            .layer_commitments
            .iter()
            .skip(1)
            .map(|commitment| {
                let beta = self.transcript.sample_field_element();
                self.transcript.append_bytes(commitment);
                beta
            })
            .collect();

        self.transcript
            .append_bytes(&proof.last_layer_value.as_bytes());

        let query_indices = (0..proof.query_decommitments.len())
            .map(|_| {
                (u64::from_be_bytes(self.transcript.sample()[..8].try_into().unwrap())
                    % self.params.domain_0_size() as u64) as usize
            })
            .collect();

        println!("  Re-created challenges by replaying transcript.");
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
    ) -> Result<(), FriValidationError> {
        let mut current_idx = query_idx;

        // Verify Merkle Paths for each layer
        for i in 0..proof.layer_commitments.len() {
            let domain_size = self.params.domain_0_size() >> i;
            let sym_idx = (current_idx + domain_size / 2) % domain_size;

            let proof_path = Proof {
                merkle_path: decommitment.layer_auth_paths[i].clone(),
            };
            if !proof_path.verify::<FriBackend>(
                &proof.layer_commitments[i],
                current_idx,
                &decommitment.layer_evaluations[i],
            ) {
                return Err(FriValidationError::InvalidMerkleProof);
            }

            let proof_path_sym = Proof {
                merkle_path: decommitment.layer_auth_paths_sym[i].clone(),
            };
            if !proof_path_sym.verify::<FriBackend>(
                &proof.layer_commitments[i],
                sym_idx,
                &decommitment.layer_evaluations_sym[i],
            ) {
                return Err(FriValidationError::InvalidMerkleProof);
            }

            println!(
                "    Layer {}: Merkle proofs valid for indices {} and {}.",
                i, current_idx, sym_idx
            );
            current_idx %= (domain_size / 2).max(1);
        }

        // Verify folding consistency from last layer to first
        let mut last_layer_eval = proof.last_layer_value.clone();
        for i in (0..proof.layer_commitments.len() - 1).rev() {
            let domain_size = self.params.domain_0_size() >> i;
            let current_query_idx = query_idx % domain_size;

            let y = &decommitment.layer_evaluations[i];
            let y_sym = &decommitment.layer_evaluations_sym[i];

            let x = generator.pow(1_u64 << i as u64).pow(current_query_idx);
            let x_inv = x.inv().unwrap();
            let two_inv = FE::from(2).inv().unwrap();

            // f_next(x^2) = (f(x) + f(-x))/2 + beta * (f(x) - f(-x))/(2x)
            let f_even = (y + y_sym) * &two_inv;
            let f_odd = (y - y_sym) * &two_inv * &x_inv;

            let expected_eval = &f_even + &betas[i] * &f_odd;

            if last_layer_eval != expected_eval {
                return Err(FriValidationError::InconsistentFolding {
                    layer: i,
                    expected: expected_eval.representative().to_hex(),
                    got: last_layer_eval.representative().to_hex(),
                });
            }
            println!("    Layer {}->{}: Folding is consistent.", i, i + 1);
            last_layer_eval = y.clone();
        }

        Ok(())
    }
}

fn run_fri_protocol() {
    // SETUP
    let poly = Polynomial::new(&[FE::from(2), -FE::from(3), FE::from(0), FE::from(1)]);
    let claimed_degree = 3;
    let params = FriParameters::new(claimed_degree, 4, 1);

    // PROVE
    let mut prover = Prover::new(poly, params.clone());
    let proof = prover.prove().unwrap();

    // VERIFY
    let mut verifier = Verifier::new(params);
    match verifier.verify(&proof) {
        Ok(_) => println!("\n✅ Proof verified successfully!"),
        Err(e) => println!("\n❌ Proof verification failed: {:?}", e),
    }
}

fn main() {
    run_fri_protocol();
}

use lambdaworks_crypto::merkle_tree::merkle::MerkleTree;
use lambdaworks_math::fft::cpu::roots_of_unity::get_powers_of_primitive_root;
use lambdaworks_math::field::traits::RootsConfig;

use crate::{FriBackend, F, FE};

/// Shared parameters for the FRI protocol, agreed upon by the Prover and Verifier.
#[derive(Debug, Clone)]
pub struct FriParameters {
    /// The initial evaluation domain (LDE).
    pub domain: Vec<FE>,
    /// How many queries the Verifier will make to check the proof.
    pub num_queries: usize,
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
        let root_order = domain_size.trailing_zeros() as u64;

        let domain =
            get_powers_of_primitive_root::<F>(root_order, domain_size, RootsConfig::Natural)
                .unwrap();

        Self {
            domain,
            num_queries,
        }
    }
}

/// Represents a single layer in the FRI protocol's commitment-folding process.
#[derive(Clone)]
pub struct FriLayer {
    /// The evaluations of the polynomial for this layer.
    pub evaluations: Vec<FE>,
    /// The Merkle tree committing to the evaluations.
    pub merkle_tree: MerkleTree<FriBackend>,
    /// The domain over which the evaluations were made.
    pub domain: Vec<FE>,
}

/// A decommitment for a single query, providing evaluations and Merkle paths for each layer.
#[derive(Debug, Clone)]
pub struct QueryDecommitment {
    /// The evaluation at the query index `q` for each layer.
    pub layer_evaluations: Vec<FE>,
    /// The Merkle authentication path for `layer_evaluations` at each layer.
    pub layer_auth_paths: Vec<Vec<[u8; 32]>>,
    /// The evaluation at the symmetric index `-q` for each layer.
    pub layer_evaluations_sym: Vec<FE>,
    /// The Merkle authentication path for `layer_evaluations_sym` at each layer.
    pub layer_auth_paths_sym: Vec<Vec<[u8; 32]>>,
}

/// The complete FRI proof sent from the Prover to the Verifier.
#[derive(Debug, Clone)]
pub struct FriProof {
    /// The Merkle root of each FRI layer.
    pub layer_commitments: Vec<[u8; 32]>,
    /// The value of the final, constant polynomial.
    pub last_layer_value: FE,
    /// The decommitments for each query.
    pub query_decommitments: Vec<QueryDecommitment>,
}

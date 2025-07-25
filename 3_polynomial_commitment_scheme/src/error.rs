use std::fmt;

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

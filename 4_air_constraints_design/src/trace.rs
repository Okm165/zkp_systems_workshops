// ================================================================================================
// PART 1: COMPUTATION & EXECUTION TRACE
// ================================================================================================
// A STARK proof starts with a computation. We represent this computation as an "execution trace".
// For this example, our computation is the Fibonacci sequence.

use crate::FE;

/// Generates a Fibonacci sequence trace of a given length.
pub fn generate_fibonacci_trace(trace_length: usize) -> Vec<FE> {
    let mut trace = vec![FE::zero(); trace_length];
    // Set the initial values, which act as our boundary conditions.
    trace[0] = FE::one();
    trace[1] = FE::one();

    for i in 2..trace_length {
        trace[i] = trace[i - 1] + trace[i - 2];
    }
    trace
}

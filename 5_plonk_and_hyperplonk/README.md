# **Chapter 5: Plonkish Arithmetization and the Shift to Multivariate Systems**

**Abstract:** This chapter introduces Plonkish arithmetization, a powerful and widely adopted alternative to the AIR framework discussed previously. We will deconstruct the core components of the Plonk protocol, including its universal gate constraint and the elegant permutation argument used to enforce wire consistency. The chapter provides a rigorous mathematical breakdown of these mechanisms and analyzes the practical efficiency of the resulting univariate PIOP. We then explore the evolution of this paradigm with HyperPlonk, a system that fundamentally shifts from univariate polynomials over multiplicative subgroups to Multi-Linear Extensions (MLEs) over the Boolean hypercube.

**Learning Objectives:** Upon completion of this chapter, you will be able to:

1.  Describe the structure of Plonkish arithmetization, including the universal gate and selector polynomials.
2.  Formally explain the permutation argument, detailing how the accumulator polynomial `Z(X)` enforces copy constraints.
3.  Analyze the prover's workload in Plonk, differentiating between the costs of NTTs and MSMs.
4.  Define the architectural shift from univariate (Plonk) to multivariate (HyperPlonk) systems, including the roles of MLEs and the Sum-check protocol.

---

## **Part 1: The Plonkish Arithmetization Paradigm**

While Chapter 4 focused on the AIR framework, characterized by uniform transition constraints well-suited for state-machine-like computations, Plonk introduces a more flexible model designed for general arithmetic circuits. It provides a structured approach to arithmetization that has become the foundation for a vast ecosystem of ZKP systems.

### **1.1 The Universal Gate Constraint**

Plonk's core innovation is a single, universal polynomial equation that can express any standard arithmetic gate. For each gate `i` in a circuit, this constraint enforces the relationship between its left (`a_i`), right (`b_i`), and output (`c_i`) wire values.

**Definition (Universal Gate Constraint).** For each gate `i ∈ {0, ..., n-1}`, the constraint is:

$$
q_{L,i} a_i + q_{R,i} b_i + q_{M,i} a_i b_i + q_{O,i} c_i + q_{C,i} = 0
$$

- The `q_{•,i}` values are **selector constants**. By enabling or disabling them, this single formula can be configured to represent:
  - **Addition (`aᵢ + bᵢ = cᵢ`):** Set `q_{L,i}=1, q_{R,i}=1, q_{O,i}=-1`.
  - **Multiplication (`aᵢ ⋅ bᵢ = cᵢ`):** Set `q_{M,i}=1, q_{O,i}=-1`.
  - **Public Constant (`aᵢ = k`):** Set `q_{L,i}=1, q_{C,i}=-k`.

### **1.2 From Gates to a Single Polynomial Identity**

To create a single, succinct proof, these `n` individual gate constraints are unified into one algebraic assertion.

1.  **Define the Domain:** We select a multiplicative subgroup `H = {ω⁰, ω¹, ..., ωⁿ⁻¹}` of a finite field `F`, where `ω` is a primitive n-th root of unity. Each point `ωⁱ` corresponds to the `i`-th gate.
2.  **Interpolate:** The wire values (`a`, `b`, `c`) and selector values (`q_L`, `q_M`, etc.) are interpolated into polynomials `a(X), b(X), c(X)` and `q_L(X), q_M(X), ...` respectively. This is practically achieved with the `O(N log N)` INTT algorithm.

**Theorem (Gate Constraint Polynomial Identity).** The `n` gate constraints are satisfied if and only if there exists a **quotient polynomial `t(X)`** such that the following identity holds for all `X`:

$$
a(X)q_L(X) + b(X)q_R(X) + a(X)b(X)q_M(X) + c(X)q_O(X) + q_C(X) = Z_H(X) \cdot t(X)
$$

Here, `Z_H(X) = X^n - 1` is the **vanishing polynomial** of the domain `H`. The existence of `t(X)` proves that the left-hand side is zero at every point in `H`.

---

## **Part 2: The Permutation Argument for Wire Consistency**

The gate constraint identity ensures local correctness but does not enforce the connections _between_ gates (i.e., that the output of one gate is the input to another). Plonk solves this "wiring problem" with its celebrated **permutation argument**.

### **2.1 The Grand Product Check**

The core idea is to prove that the set of all `3n` wire values is a permutation of itself according to the circuit's wiring diagram.

1.  **Define Permutation Polynomials:** The wiring is encoded into pre-computed **permutation polynomials** `S_{σ,1}(X), S_{σ,2}(X), S_{σ,3}(X)`. For a wire `(j, i)` (e.g., wire `aᵢ`), `S_{σ,j}(ωⁱ)` evaluates to the unique global identity of the wire it is connected to.
2.  **Randomized Binding:** The verifier provides random challenges `β` and `γ`. These are used to "bind" each wire's value to its unique identity, preventing fraudulent swaps.
3.  **The Accumulator `Z(X)`:** Instead of checking a massive product directly, Plonk uses a recursive **accumulator polynomial `Z(X)`**. This polynomial is constrained by two identities that must hold over `H`:
    - **Start Condition:** The accumulator must start at 1. `L_1(X)` is the first Lagrange polynomial.
      $$
      L_1(X) \cdot (Z(X) - 1) = 0
      $$
    - **Recursive Step:** The core logic. `Z(Xω)` must equal `Z(X)` multiplied by the ratio of the randomized "original" wire values to the "permuted" wire values. This is rearranged to form a polynomial constraint:
      $$
      Z(X\omega) \prod_{j=1}^{3}(w_j(X) + \beta S_{\sigma,j}(X) + \gamma) - Z(X) \prod_{j=1}^{3}(w_j(X) + \beta S_{id,j}(X) + \gamma) = 0
      $$
      The prover must demonstrate that this polynomial, along with the start condition polynomial, is also divisible by `Z_H(X)`.

---

## **Part 3: The Shift to Multivariate Systems - HyperPlonk**

HyperPlonk was designed to address a perceived bottleneck in Plonk: the `O(N log N)` complexity of the FFT. It achieves this by fundamentally changing the underlying mathematical representation of the proof.

### **3.1 From Univariate to Multi-Linear Extensions (MLEs)**

The core architectural shift is from univariate polynomials over multiplicative subgroups to **Multi-Linear Extensions (MLEs)** over the **Boolean hypercube**.

**Definition (Multi-Linear Extension, MLE).** For a function `f: \{0,1\}^v \to \mathbb{F}`, its unique MLE `f̃` is the `v`-variate polynomial of degree 1 in each variable that agrees with `f` on all `2^v` points of the Boolean hypercube `{0,1\}^v`.

$$
\tilde{f}(X_1, \dots, X_v) = \sum_{w \in \{0,1\}^v} f(w) \cdot \prod_{i=1}^{v} (X_i w_i + (1-X_i)(1-w_i))
$$

### **3.2 From Quotients to the Sum-check Protocol**

This change of representation makes the quotient-based checks of Plonk inapplicable. HyperPlonk replaces them with the **Sum-check protocol**, an interactive protocol that allows a prover to convince a verifier of the value of a sum `S = ∑_{x∈\{0,1\}^v} g(x)` with `O(log N)` verifier complexity.

- **Gate Constraints via ZeroCheck:** To prove a constraint `C(X) = 0` holds over the hypercube, the verifier sends a random point `r`, and the prover uses Sum-check to prove:
  $$
  \sum_{x \in \{0,1\}^v} C(x) \cdot \text{eq}(x, r) = 0
  $$
  This reduces the check to `C(r)=0`, which implies `C(X)` is the zero polynomial with high probability.
- **Permutation Constraints via Log-Derivative:** The permutation argument is transformed into a sum that must equal zero, which is then verified with Sum-check. For a random challenge `γ`:
  $$
  \sum_{x \in \{0,1\}^v} \left( \frac{1}{\tilde{f}(x) + \gamma} - \frac{1}{\tilde{g}(x) + \gamma} \right) = 0
  $$

---

## **Part 4: Comparative Analysis: Plonk vs. HyperPlonk**

This section provides a rigorous analysis of the practical trade-offs between the two systems.

| Feature                     | Plonk (Univariate)                                                                                   | HyperPlonk (Multivariate)                                                                      |
| :-------------------------- | :--------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------- |
| **Prover Complexity**       | `O(N log N)`                                                                                         | `O(N)`                                                                                         |
| **Practical Bottleneck**    | **Compute-bound.** Dominated by `O(N)` expensive MSMs, not the `O(N log N)` NTT.                     | **Memory-bound.** Dominated by the serial, high-memory-traffic Sum-check protocol.             |
| **Hardware Friendliness**   | **High.** NTTs are highly parallelizable and well-optimized in hardware.                             | **Low.** The serial nature of Sum-check is difficult to accelerate effectively.                |
| **Verifier Complexity**     | **`O(1)` (Constant).** Requires 2 pairings.                                                          | **`O(log N)` (Logarithmic).** Requires `log N` pairings.                                       |
| **On-Chain Viability**      | **Excellent.** Low, predictable gas cost.                                                            | **Poor.** Prohibitively high gas cost for large circuits.                                      |
| **Proof Size**              | **`O(1)` (Constant).** ~0.5-1 KB.                                                                    | **`O(log N)` (Logarithmic).** ~5-10 KB.                                                        |
| **Custom Gate Flexibility** | **Limited.** High-degree gates increase the degree of the quotient polynomial, raising prover costs. | **High.** The Sum-check protocol is tolerant of high-degree constraints with minimal overhead. |

---

**Author:** [Okm165](https://github.com/Okm165) | [@bartolomeo_diaz](https://x.com/bartolomeo_diaz)

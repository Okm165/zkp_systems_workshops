# Zero-Knowledge Proof systems: A Deep Dive from Foundations to Frontiers

Welcome to a comprehensive course on the theory and practice of modern Zero-Knowledge Proof (ZKP) systems. This repository contains a series of lectures designed to deconstruct the complex world of ZKPs, starting from the mathematical first principles and building up to the architectural design of state-of-the-art systems.

---

## Course Lectures

The course is structured as a series of chapters, each building on the last.

- ### [Chapter 0: The Architecture of a Zero-Knowledge Proof System](./0_zkp_architecture/README.md)

  - **Description:** This chapter introduces the "big picture" of a modern ZKP system. We define the axiomatic guarantees of a proof system (Completeness, Soundness, Zero-Knowledge) and present the powerful architectural framework that separates the logic layer (PIOPs) from the cryptographic layer (PCSs).

- ### [Chapter 1: The Mathematical Toolkit for Verifiable Computation](./1_mathematical_toolkit/README.md)

  - **Description:** Here, we establish the three mathematical pillars of ZKPs: Finite Fields for perfect arithmetic, Polynomials for encoding computational logic (Arithmetization), and Cryptographic Hash Functions for building transparent and non-interactive systems.

- ### [Chapter 2: The Fast Fourier Transform for Polynomial Multiplication](./2_fast_polynomial_arithmetic/README.md)

  - **Description:** This chapter introduces the Fast Fourier Transform (FFT) as a high-performance algorithm for polynomial multiplication. We analyze the $`\Theta(n^2)`$ bottleneck of naive convolution and demonstrate how the FFT, by using the complex roots of unity, provides a $`\Theta(n \log n)`$ method for converting polynomials to and from the point-value representation, where multiplication is a linear-time operation. We derive the algorithm's structure, including the butterfly operation and bit-reversal permutation, and show how the same logic applies to the inverse transform.

- ### [Chapter 3: Foundations of Polynomial Commitment Schemes with a Focus on FRI](./3_polynomial_commitment_scheme/README.md)

  - **Description:** This chapter introduces Polynomial Commitment Schemes (PCS) as a foundational cryptographic primitive essential for modern Zero-Knowledge Proof (ZKP) systems. We will define a PCS and its core properties—binding, succinctness, and evaluation proofs—and explore its role in translating computational integrity claims into verifiable algebraic statements

- ### [Chapter 4: Algebraic Intermediate Representation (AIR) and Constraint Design](./4_air_constraints_design/README.md)
  - **Description:** This chapter introduces the Algebraic Intermediate Representation (AIR), the framework for converting computational claims into polynomial constraints. We detail the process of arithmetization by designing boundary and transition constraints, unifying them into a composition polynomial, and demonstrating how out-of-domain checks and the DEEP composition polynomial work together to ensure the integrity of a STARK proof.

---

## About & Contributions

Contributions, feedback, and questions are highly welcome. Please feel free to open an issue or a pull request to help improve the material for everyone.

**Author:** [Okm165](https://github.com/Okm165) | [@bartolomeo_diaz](https://x.com/bartolomeo_diaz)

## License

This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/).

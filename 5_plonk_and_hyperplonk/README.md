# **Chapter 5: Plonkish Arithmetization and the Shift to Multivariate Systems**

**Abstract:** This chapter introduces Plonkish arithmetization, a powerful and widely adopted alternative to the AIR framework discussed previously. We will deconstruct the core components of the Plonk protocol, including its universal gate constraint and the elegant permutation argument used to enforce wire consistency. The chapter provides a rigorous mathematical breakdown of these mechanisms and analyzes the practical efficiency of the resulting univariate PIOP. We then explore the evolution of this paradigm with HyperPlonk, a system that fundamentally shifts from univariate polynomials over multiplicative subgroups to Multi-Linear Extensions (MLEs) over the Boolean hypercube. This transition, motivated by the pursuit of a linear-time prover, replaces the FFT and quotient checks of Plonk with the Sum-check protocol.

**Learning Objectives:** Upon completion of this chapter, you will be able to:

1.  Describe the structure of Plonkish arithmetization, including the universal gate and selector polynomials.
2.  Formally explain the permutation argument, detailing how the accumulator polynomial $Z(X)$ enforces copy constraints.
3.  Define the architectural shift from univariate (Plonk) to multivariate (HyperPlonk) systems, including the roles of MLEs and the Sum-check protocol.

---

## **Part 1: The Plonkish Arithmetization Paradigm**

While Chapter 4 focused on the AIR framework, characterized by uniform transition constraints well-suited for state-machine-like computations, Plonk introduces a more flexible model designed for general arithmetic circuits. It provides a structured approach to arithmetization that has become the foundation for a vast ecosystem of ZKP systems.

### **1.1 The Universal Gate Constraint**

Plonk's core innovation is a single, universal polynomial equation that can express any standard arithmetic gate. For each gate $i$ in a circuit, this constraint enforces the relationship between its left ($a_i$), right ($b_i$), and output ($c_i$) wire values.

**Definition (Universal Gate Constraint).** For each gate $i \in \{0, \dots, n-1\}$, the constraint is:

$$
\begin{align*}
q_{L,i} a_{i} + q_{R,i} b_{i} + q_{M,i} a_{i} b_{i} + q_{O,i} c_{i} + q_{C,i} = 0
\end{align*}
$$

- The $q_{\bullet,i}$ values are **selector constants**. By enabling or disabling them, this single formula can be configured to represent different operations.

### **1.2 From Gates to a Single Polynomial Identity**

To create a single, succinct proof, these $n$ individual gate constraints are unified into one algebraic assertion. This is done by interpolating the wire and selector values into polynomials over a multiplicative subgroup $H = \{\omega^{0}, \dots, \omega^{n-1}\}$ of a finite field $\mathbb{F}$.

**Theorem (Gate Constraint Polynomial Identity).** The $n$ gate constraints are satisfied if and only if there exists a **quotient polynomial $t(X)$** such that the following identity holds:

$$
\begin{align*}
a(X)q_L(X) + b(X)q_R(X) + a(X)b(X)q_M(X) + c(X)q_O(X) + q_C(X) = Z_H(X) \cdot t(X)
\end{align*}
$$

Here, $Z_H(X) = X^n - 1$ is the **vanishing polynomial** of the domain $H$.

---

## **Part 2: The Permutation Argument for Wire Consistency (Expanded)**

The gate constraint identity from Part 1 is a necessary but insufficient condition for a valid proof. It confirms local correctness but fails to enforce the connections _between_ gates. This section provides a detailed, step-by-step deconstruction of Plonk's solution: the permutation argument.

### **2.1 The High-Level Idea: A Permutation Check**

The circuit's wiring diagram defines a specific **permutation**, $\sigma$, on the _positions_ of all $3n$ wires. If all copy constraints are satisfied, then reordering the _values_ according to this permutation $\sigma$ should result in an identical collection of values.

### **2.2 The Tool: Checking Multiset Equality with a Grand Product**

To check if two multisets $\{f_i\}$ and $\{g_i\}$ are equal, we check a randomized product identity. For a random challenge $\gamma$, the sets are equal if and only if:

$$
\begin{align*}
\prod_{i=1}^{k} (f_i + \gamma) = \prod_{i=1}^{k} (g_i + \gamma)
\end{align*}
$$

This is the **grand product check**. To enforce a specific permutation and prevent fraudulent swaps, we enhance this check with another random challenge, $\beta$, which binds each wire's value to its unique identity.

### **2.3 The Formal Grand Product Identity of Plonk**

The prover must demonstrate that the multiset of randomized wire values is invariant when we swap the wire identities ($S_{id,j}$) with the permuted wire identities ($S_{\sigma,j}$).

$$
\begin{align*}
\prod_{i=0}^{n-1} \prod_{j=1}^{3} (w_j(\omega^i) + \beta S_{id,j}(\omega^i) + \gamma) = \prod_{i=0}^{n-1} \prod_{j=1}^{3} (w_j(\omega^i) + \beta S_{\sigma,j}(\omega^i) + \gamma)
\end{align*}
$$

### **2.4 The Accumulator $Z(X)$: Efficient Verification**

Computing the two giant products is inefficient. Instead, we check that their ratio is 1 by building the ratio step-by-step using a recursive **accumulator polynomial, $Z(X)$**.

1.  **Recursive Relation:** We define $Z(X)$ such that it starts at 1 ($Z(\omega^0)=1$) and is updated at each step:
$$
Z(\omega^{i+1}) = Z(\omega^i) \cdot \frac{\prod_{j=1}^{3}(w_j(\omega^i) + \beta S_{\sigma,j}(\omega^i) + \gamma)}{\prod_{j=1}^{3}(w_j(\omega^i) + \beta {id,j}(\omega^i) + \gamma)}
$$
2.  **Final Condition:** If the grand product identity holds, the accumulator must return to its starting value: $Z(\omega^n) = 1$.

### **2.5 Finalizing the Polynomial Constraints**

These conditions on $Z(X)$ are converted into polynomial constraints that must be divisible by $Z_H(X)$:

- **Start Constraint:** $L_1(X) \cdot (Z(X) - 1) = 0$, where $L_1(X)$ is the first Lagrange polynomial.
- **Recursive Step Constraint:**
  ```math
    \begin{align*}
    Z(X\omega) \prod_{j=1}^{3}(w_j(X) + \beta S_{id,j}(X) + \gamma) - Z(X) \prod_{j=1}^{3}(w_j(X) + \beta S_{\sigma,j}(X) + \gamma) = 0
    \end{align*}
  ```
  These two new constraints are bundled with the gate constraints into the final quotient polynomial t(X).

---

## **Part 3: The Shift to Multivariate Systems - HyperPlonk (Expanded)**

HyperPlonk was designed to attack the theoretical $O(N \log N)$ complexity of the FFT. To achieve this, it reimagines the two pillars of Plonk's arithmetization: its **domain** and its **polynomial representation**.

### **3.1 Step 1: A New Domain - The Boolean Hypercube**

To eliminate the FFT, HyperPlonk chooses the **Boolean hypercube**, $\{0,1\}^v$, where $N=2^v$. The execution trace is re-indexed so that the value at gate $i$ corresponds to the value at the vertex represented by the binary form of $i$.

### **3.2 Step 2: A New Representation - Multi-Linear Extensions (MLEs)**

A univariate polynomial is no longer a natural fit. We need a polynomial in $v$ variables. HyperPlonk uses a specific, highly structured type: the **Multi-Linear Extension (MLE)**.

**Definition (Multi-Linear Extension, MLE).** For a function $f: \{0,1\}^v \to \mathbb{F}$, its unique MLE $\tilde{f}$ is the $v$-variate polynomial of degree 1 in each variable that agrees with $f$ on all $2^v$ points of the Boolean hypercube $\{0,1\}^v$.

$$
\begin{align*}
\tilde{f}(X_1, \dots, X_v) = \sum_{w \in \{0,1\}^v} f(w) \cdot \prod_{i=1}^{v} (X_i w_i + (1-X_i)(1-w_i))
\end{align*}
$$

The product term, $\text{eq}(X, w)$, acts as a selector that is 1 only when the input $X$ equals the corner $w$.

### **3.3 Step 3: A New Tool - The Sum-check Protocol**

This new representation makes Plonk's quotient checks inapplicable. HyperPlonk introduces the **Sum-check protocol**, an efficient interactive protocol to prove the value of a sum $S = \sum_{x \in \{0,1\}^v} g(x)$ with only $O(\log N)$ verifier work.

### **3.4 Applying the Tool: Constraint Verification**

All of HyperPlonk's constraint checks are reduced to problems that the Sum-check protocol can solve.

#### **Gate Constraints via ZeroCheck**

To prove a constraint $C(X)$ is zero on all $2^v$ points, the verifier sends a random point $r$, and the prover uses Sum-check to prove:

$$
\begin{align*}
\sum_{x \in \{0,1\}^v} C(x) \cdot \text{eq}(x, r) = 0
\end{align*}
$$

Since $\text{eq}(x, r)$ isolates the term $C(r)$, this is equivalent to proving $C(r)=0$.

#### **Permutation Constraints via Log-Derivative**

To adapt Plonk's grand product check into a sum, the prover uses Sum-check to prove the following sum is zero for a random $\gamma$:

$$
\begin{align*}
\sum_{x \in \{0,1\}^v} \left( \frac{1}{\tilde{f}(x) + \gamma} - \frac{1}{\tilde{g}(x) + \gamma} \right) = 0
\end{align*}
$$

Here, $\tilde{f}$ represents the original randomized wire values and $\tilde{g}$ the permuted ones.

---

**Author:** [Okm165](https://github.com/Okm165) | [@bartolomeo_diaz](https://x.com/bartolomeo_diaz)

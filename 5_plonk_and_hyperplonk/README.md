# **Chapter 5: Plonkish Arithmetization and the Shift to Multivariate Systems**

**Abstract:** This chapter introduces Plonkish arithmetization, a powerful and widely adopted alternative to the AIR framework discussed previously. We will deconstruct the core components of the Plonk protocol, including its universal gate constraint and the elegant permutation argument used to enforce wire consistency. The chapter provides a mathematical breakdown of these mechanisms and analyzes the practical efficiency of the resulting univariate PIOP. We then explore the evolution of this paradigm with HyperPlonk, a system that fundamentally shifts from univariate polynomials over multiplicative subgroups to Multi-Linear Extensions (MLEs) over the Boolean hypercube. This transition, motivated by the pursuit of a linear-time prover, replaces the FFT and quotient checks of Plonk with the Sum-check protocol.

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

The gate constraint identity from Part 1 is a necessary but insufficient condition for a valid proof. It confirms that each gate performs its operation correctly in isolation (e.g., this gate adds, that gate multiplies). However, it fails to enforce the connections _between_ the gates. This section provides a detailed, step-by-step deconstruction of Plonk's solution: the permutation argument.

### **2.1 The High-Level Idea: A Permutation Check**

Imagine we list all $3n$ wire values in our circuit in one large collection. The circuit's wiring diagram defines a specific **permutation**, $\sigma$, on the _positions_ of these wires. For instance, if the wire at position $c_0$ (output of gate 0) is connected to the wire at position $a_1$ (left input of gate 1), the permutation maps these positions to each other.

If all copy constraints are satisfied (i.e., $c_0 = a_1$), then if we were to reorder the _values_ according to this permutation $\sigma$, the resulting collection of values would be identical to the original. The permutation argument's goal is to prove this property efficiently and in zero-knowledge.

### **2.2 The First Tool: Checking Multiset Equality with a Grand Product**

How can we prove that two multisets of numbers, $F = \{f_1, \dots, f_k\}$ and $G = \{g_1, \dots, g_k\}$, are permutations of each other? A simple sum or product check is insufficient. The cryptographic solution is to use a randomized product.

For a random challenge $\gamma$ provided by the verifier, the sets are equal if and only if the following identity holds with overwhelming probability:

$$
\begin{align*}
\prod_{i=1}^{k} (f_i + \gamma) = \prod_{i=1}^{k} (g_i + \gamma)
\end{align*}
$$

By the Schwartz-Zippel lemma, this polynomial identity in $\gamma$ implies that the underlying multisets are identical. This is the **grand product check**.

### **2.3 Refining the Argument: Enforcing a Specific Permutation**

The simple grand product check only confirms that the set of wire values is the same before and after the permutation. It doesn't enforce the specific connections. Furthermore, it doesn't distinguish between a value $v$ appearing in wire column $a$ versus the same value $v$ appearing in column $b$.

To solve this, we enhance the check with more randomness. The verifier provides a second random challenge, $\beta$. We then give each wire position a unique identity and bind it to its value.

1.  **Unique Wire Identities:** Each of the $3n$ wire positions is given a unique identifier. In Plonk, these are represented by pre-computed **identity polynomials**, $S_{id,1}(X), S_{id,2}(X), S_{id,3}(X)$, where $S_{id,j}(\omega^i)$ gives the unique ID for the wire at gate $i$, column $j$.
2.  **Permutation Polynomials:** The wiring permutation $\sigma$ is also encoded into polynomials, $S_{\sigma,1}(X), S_{\sigma,2}(X), S_{\sigma,3}(X)$. Here, $S_{\sigma,j}(\omega^i)$ gives the ID of the wire that $(j, i)$ is connected to.
3.  **Randomized Binding:** The value at each wire is now combined with its identity using $\beta$ and $\gamma$. For a wire $(j, i) with value $w_j(\omega^i)$ and identity $S_{id,j}(\omega^i)$, its randomized value is $w_j(\omega^i) + \beta \cdot S_{id,j}(\omega^i) + \gamma$.

Now, the core assertion is that the multiset of randomized values is invariant when we swap the wire identities with the permuted wire identities. This leads to the **final grand product identity**:

$$
\begin{align*}
\prod_{i=0}^{n-1} \prod_{j=1}^{3} (w_j(\omega^i) + \beta S_{id,j}(\omega^i) + \gamma) = \prod_{i=0}^{n-1} \prod_{j=1}^{3} (w_j(\omega^i) + \beta S_{\sigma,j}(\omega^i) + \gamma)
\end{align*}
$$

This identity holds if and only if all copy constraints are satisfied. If even one wire has an incorrect value (e.g., $c_0 \neq a_1$), then the left and right products will differ, and the equality will fail with high probability due to the randomness of $\beta$ and $\gamma$.

### **2.5 The Accumulator $Z(X)$: Making the Check Efficient**

Computing two giant products of $3n$ terms is inefficient. Plonk uses a clever trick to check this identity recursively. The equality of the two grand products is equivalent to their ratio being 1. We can build this ratio step-by-step using an **accumulator polynomial, $Z(X)$**.

1.  **Define a Recursive Relation:** We define $Z(X)$ such that it starts at 1 and, at each step `i`, it is multiplied by the ratio of that step's contribution to the two products.
    _ **Start:** $Z(\omega^0) = 1$
    _ **Recurrence:** For $i = 0, \dots, n-1$:
    $$
    \begin{align*}
    Z(\omega^{i+1}) = Z(\omega^i) \cdot \frac{\prod_{j=1}^{3}(w_j(\omega^i) + \beta S_{\sigma,j}(\omega^i) + \gamma)}{\prod_{j=1}^{3}(w_j(\omega^i) + \beta S_{id,j}(\omega^i) + \gamma)}
    \end{align*}
    $$
2.  **Final Condition:** If the grand product identity holds, the product of all the numerators will equal the product of all the denominators over the full cycle. Therefore, the accumulator must return to its starting value: $Z(\omega^n) = Z(\omega^0) = 1$.

### **2.6 Finalizing the Polynomial Constraints: From Recurrence to Polynomials**

This is the crucial conceptual leap. The conditions on $Z(X)$ are defined at discrete points ($\omega^i$, $\omega^{i+1}$), but our proof system requires constraints that are **polynomials**. We must convert these discrete checks into polynomial identities that are zero everywhere on our domain $H$.

#### **Constraint 1: The Starting Value**

- **The Condition:** We need $Z(\omega^0) = 1$.
- **The Polynomial Equivalent:** This is equivalent to saying that the polynomial $Z(X) - 1$ must have a root at $X = \omega^0$.
- **Enforcing the Constraint:** To enforce this without affecting other points, we multiply by the first Lagrange basis polynomial, $L_1(X)$, which is `1` at $\omega^0$ and `0` at all other points in $H$. This gives us our first constraint polynomial, which must be zero on all of $H$:
  $$
  \begin{align*}
  L_1(X) \cdot (Z(X) - 1) = 0, \quad \forall X \in H
  \end{align*}
  $$

#### **Constraint 2: The Recursive Step**

- **The Condition:** For each $i \in \{0, \dots, n-2\}$, the recurrence relation must hold. Let's first rearrange the recurrence algebraically to eliminate the fraction:
  $Z(\omega^{i+1}) \cdot (\text{denominator at } i) = Z(\omega^i) \cdot (\text{numerator at } i)$
- **The Polynomial Equivalent:** Now we "lift" this discrete relation into the world of polynomials.
  - The value at the "current" step $i$, like $Z(\omega^i)$, is simply the polynomial $Z(X)$ evaluated at $X=\omega^i$.
  - The value at the "next" step, $Z(\omega^{i+1})$, can be represented by evaluating the same polynomial $Z(X)$ at the "next" point in the domain, which is $X \cdot \omega$. So, $Z(\omega^{i+1})$ is represented by $Z(X\omega)$.
- **Enforcing the Constraint:** By substituting the polynomials for their evaluations, we get a single polynomial equation that must hold for all $X \in H$ (except the very last point, which is handled by a small adjustment not detailed here). This is our second constraint polynomial:
  $$
  \begin{align*}
  Z(X\omega) \prod_{j=1}^{3}(w_j(X) + \beta S_{id,j}(X) + \gamma) - Z(X) \prod_{j=1}^{3}(w_j(X) + \beta S_{\sigma,j}(X) + \gamma) = 0, \quad \forall X \in H
  \end{align*}
  $$
        *Notice the terms are swapped compared to the fraction: the denominator from the recurrence is now multiplied by $Z(X\omega)$, and the numerator is multiplied by $Z(X)$.*

These two new polynomial constraints, which must be divisible by $Z_H(X)$, now perfectly capture the logic of the accumulator. They are bundled with the gate constraints using random challenges into the final, single quotient polynomial $t(X)$. This elegant construction transforms a complex, global graph property (circuit wiring) into a set of local, algebraic polynomial constraints.

---

## **Part 3: The Shift to Multivariate Systems - HyperPlonk**

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

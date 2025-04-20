# 🧠 Libint-style Integral Evaluation Documentation

This file documents the implementation of **1-electron** and **2-electron** integrals using Obara-Saika methods in Hartree-Fock theory.

---

## 📘 `1e.cpp` – One-Electron Integrals

```cpp
// Core integrals implemented:
//  - Overlap (S)
//  - Kinetic (T)
//  - Nuclear attraction (V)
//  - Obara-Saika recursion helper: E(...)
//  - Boys function recursion helper: R(...)
//  - Gaussian product center: P(...)
```

### 🔹 `double Sp(...)`
Evaluates **primitive overlap** using the Obara-Saika formula.

\[
S_{ab} = E(l_a, l_b, 0, \alpha_a, \alpha_b, A_x - B_x) \cdot E(...) \cdot E(...) \cdot \left( \frac{\pi}{\alpha_a + \alpha_b} \right)^{3/2}
\]

### 🔹 `double S(GF g1, GF g2)`
Computes **contracted overlap** between `g1` and `g2`:

\[
S = \sum_{i,j} N_i N_j d_i d_j S_{ij}
\]

### 🔹 `double Tp(...)`
Computes **primitive kinetic energy integral** via Obara-Saika recursion:

\[
T = -\frac{1}{2} \nabla^2 \langle a | b \rangle
\]

### 🔹 `double T(GF g1, GF g2)`
Computes **contracted kinetic energy**:

\[
T = \sum_{i,j} N_i N_j d_i d_j T_{ij}
\]

### 🔹 `double R(...)`
Implements **Obara-Saika recursion** for Boys-integral-based Hermite expansions.

### 🔹 `double Vp(...)`
Evaluates **primitive nuclear attraction integral**.

### 🔹 `double V(GF g1, GF g2, std::vector<double> xyzN)`
Computes **contracted nuclear attraction**.

---

## 📘 `2e.cpp` – Two-Electron Integrals

```cpp
// Implements:
//  - 2-electron (electron repulsion) integrals ⟨ab|cd⟩
//  - Symmetry-aware caching
```

### 🔹 `ERI::isEqual()`
Tests **8-fold permutational symmetry**.

### 🔹 `double Gp(...)`
Computes **primitive ERI** using Obara-Saika:

\[
\langle ab|cd \rangle = \sum_{t,u,v,\tau,\nu,\phi} E(...) \cdot R(0, t+\tau, u+\nu, v+\phi)
\]

### 🔹 `double G(...)`
Contracted **electron repulsion** integral:

\[
G = \sum_{ijkl} N_i N_j N_k N_l \cdot d_i d_j d_k d_l \cdot G_{ijkl}
\]

### 🔹 `ERIs(...)`
Returns a 4D tensor of ⟨ab|cd⟩ integrals with symmetry.

---

## ✅ Summary Table

| Function         | Purpose                        | Equation Type      |
|------------------|--------------------------------|--------------------|
| `Sp` / `S`       | Overlap integral               | \(\langle a | b \rangle\) |
| `Tp` / `T`       | Kinetic energy                 | \(-\frac{1}{2} \langle a | \nabla^2 | b \rangle\) |
| `Vp` / `V`       | Nuclear attraction             | \(\langle a | \frac{1}{r} | b \rangle\) |
| `Gp` / `G`       | Electron-electron repulsion    | \(\langle ab | cd \rangle\) |
| `R`              | Boys-function recursion        | \(\text{Boys}(n, x)\) |
| `E`              | Hermite Gaussian expansion     | Obara-Saika term   |
| `P`              | Gaussian product center        | \(\vec{P} = \frac{\alpha_1 \vec{A} + \alpha_2 \vec{B}}{\alpha_1 + \alpha_2}\) |
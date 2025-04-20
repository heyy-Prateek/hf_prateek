
# 🧮 LibMath Documentation

This library contains mathematical helper utilities used for Hartree-Fock quantum chemistry. It is divided into:

- `linalg.hpp` / `linalg.cpp`: Matrix operations (linear algebra)
- `miscmath.hpp` / `miscmath.cpp`: Special functions (factorials, Boys function)

---

## 📘 `linalg.hpp` + `linalg.cpp` – Matrix Utilities

### `Matrix` Class

Represents a 2D matrix with support for:

- Dynamic allocation/deallocation
- Operator overloads (`+`, `-`, `*`, assignment)
- Transposition and indexing
- Diagonalization using LAPACK routines

### 🔧 Key Methods

- `Matrix(int r, int c, bool sym = false)`: Initialize matrix
- `~Matrix()`: Destructor
- `Matrix::operator*`, `+`, `-`: Element-wise and matrix-matrix operations
- `transpose(Matrix A)`: Return transpose
- `I(int n, int n)`: Identity matrix
- `zero(int r, int c)`: Zero matrix
- `Tr(Matrix A)`: Trace
- `diagonalize(Matrix A)`: Diagonalizes using LAPACK's `dsyev_`
- `m_sqrt(Matrix A)`: Matrix square root
- `m_inv_sqrt(Matrix A)`: Inverse square root
- `sym_linear_solve(Matrix A, Matrix B, int* icd)`: Solve `AX = B` using LAPACK's `dsysv_`

---

## 📘 `miscmath.hpp` + `miscmath.cpp` – Special Math

### ✅ Factorials and Double Factorials

```cpp
double long fact(double long n)
```

- Computes regular factorial `n!`

```cpp
double long dfact(double long n)
```

- Computes double factorial: `n!!` (used in Gaussian normalization and Boys function)

---

### 🔁 Dot Product

```cpp
double dot(std::vector<double> a, std::vector<double> b)
```

- Computes `a ⋅ b`

---

### 🧠 Boys Function

```cpp
double boys(int n, double x)
```

Computes the Boys function:
\[
F_n(x) = \int_0^1 t^{2n} e^{-x t^2} dt
\]

Efficient implementation based on Taylor expansion or asymptotic expansion depending on the value of `x`.

---

## 🧪 LAPACK Integration

Functions like `diagonalize()` and `sym_linear_solve()` rely on LAPACK routines:
- `dsyev_` for symmetric eigensystems
- `dsysv_` for solving linear systems

Ensure LAPACK is linked during compilation.

---

## ✅ Summary Table

| Function | Description |
|----------|-------------|
| `Matrix` | General matrix operations |
| `diagonalize` | Returns eigenvalues/eigenvectors |
| `m_sqrt` | Computes matrix square root |
| `boys` | Boys function for integrals |
| `fact` | Factorial |
| `dfact` | Double factorial |

---

Made for Hartree-Fock SCF but reusable across scientific computing!


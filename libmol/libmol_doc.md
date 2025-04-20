
# 📦 libmol – Molecule & Basis Set Handler

This file handles:
- Reading molecule files (`.inp`)
- Computing number of electrons and spin configuration
- Assigning basis sets to atoms (`STO-3G` currently hardcoded)
- Generating `GF` (Gaussian Function) objects for each atom

---

## 🔹 `Molecule::Molecule(std::string file, std::string bfs)`

### Reads input and initializes:

| Property     | Description                           |
|--------------|---------------------------------------|
| `Natoms`     | Number of atoms                       |
| `charge`     | Total molecular charge                |
| `multiplicity` | Spin multiplicity (e.g., 1 = singlet, 2 = doublet) |
| `Zvals`      | Atomic numbers for all atoms          |
| `xyz`        | Cartesian coordinates of atoms        |
| `basis`      | Basis set name (e.g., `STO-3G`)       |
| `Nelec`      | Total number of electrons             |
| `NUPDOWN`    | `multiplicity - 1` (used for UHF)     |

### 🧪 Important logic

- Computes number of electrons from atomic numbers & charge
- Calculates number of spin-up vs spin-down electrons via `NUPDOWN`
- Initializes atomic orbitals using `AOfunctions(...)`

---

## 🔸 `std::vector<GF> AOfunctions(...)`

Returns all Gaussian functions (contracted GTOs) for a given element and position.

### Input Parameters:
- `bfs` → Basis set name (`"STO-3G"`)
- `Zval` → Atomic number
- `pos` → 3D position (center of basis functions)

### Currently supported: `STO-3G` only

Each `Zval` creates appropriate `GF` instances based on hardcoded STO-3G coefficients.

---

### ✅ Example: Carbon (Z=6)

- 1s orbital
- 2s and 2p orbitals: `2px`, `2py`, `2pz`

All use exponents `α`, contraction coefficients `d`, and angular momentum (L) like:

```cpp
std::vector<double> a2sp = {2.941249355, 0.6834830964, 0.2222899159};
std::vector<double> d2s  = {-0.09996722919, 0.3995128261, 0.7001154689};
std::vector<double> d2p  = {0.1559162750, 0.6076837186, 0.3919573931};
```

These are passed to:

```cpp
GF C2s(a2sp, d2s, pos, Ls);   // s orbital
GF C2px(a2sp, d2p, pos, Lpx); // px orbital
...
```

---

## 🧠 How to Add More Basis Sets (Hardcoded)

You can **extend the current switch-case logic** in `AOfunctions(...)` like this:

```cpp
if (bfs == "STO-3G") {
    // already implemented
}
else if (bfs == "6-31G") {
    if (Zval == 1) {
        std::vector<double> a1s = {...};
        std::vector<double> d1s = {...};
        orbitals.push_back(GF(a1s, d1s, pos, Ls));
    }
    // repeat for other elements like C, N, O
}
else {
    std::cerr << "Unsupported basis set: " << bfs << std::endl;
    exit(1);
}
```

📌 Store your `6-31G`, `cc-pVDZ`, or `def2-SVP` parameters in header files or separate `*.cpp` modules if it gets too long.

---

## ⚠️ Future Improvements

- 📁 Use external basis set libraries like [Basis Set Exchange](https://www.basissetexchange.org/)
- 📜 Parse basis sets from JSON or BSE-style plain text
- ⚙️ Switch to data-driven design (no more hardcoding)
- ✅ Add support for d- and f-functions (`l + m + n > 1`)

---

## 🔚 Summary

| Component        | Role                                              |
|------------------|---------------------------------------------------|
| `Molecule`       | Initializes a molecule from an input `.inp` file  |
| `AOfunctions`    | Maps basis set + atomic number to `GF` orbitals   |
| `STO-3G` Support | Fully hardcoded for H, He, C, N, O                |
| Extendable?      | ✅ Yes, via `if (bfs == "...")` blocks             |

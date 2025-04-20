
# 🧪 Hartree–Fock Quantum Chemistry Solver By Prateek Saxena(2021B2A12756P)_BITS_Pilani

This repository contains a modular, header-based implementation of the **Hartree-Fock (HF)** method in C++, including:
- Support for **RHF** and **UHF** with DIIS and fixed-point SCF
- Mulliken and Löwdin population analysis
- Gaussian integrals via **Libint** and hardcoded STO-3G basis sets
- Modular library structure: `libmol`, `libmath`, `libgf`, `libscf`, and integral routines

---

## Contact
Email: f20212756@pilani.bits-pilani.ac.in
Contact No: 9458827686
Address: BITS Pilani_Pilani_Rajasthan_333031_India

---

## 📁 Directory Structure

```
├── main/                   # Contains `hartree-fock.cpp` and `run.sh`
├── libmol/                # Molecule class and basis set loader
├── libmath/               # Linear algebra utilities and math helpers
├── libgf/                 # Gaussian function definition and normalization
├── libscf/                # SCF driver for RHF/UHF + DIIS/FPI
├── libint/                # One- and two-electron integrals using libint
└── README.md              # This file
```

---

## 🧪 Input File Format

```
<Number of atoms> <charge> <multiplicity>
<Z1> <x1> <y1> <z1>
<Z2> <x2> <y2> <z2>
...
<Basis Set Name>
```

Example for water (H₂O):
```
3 0 1
8  0.0000 0.0000 0.0000
1  0.7580 0.5860 0.0000
1 -0.7580 0.5860 0.0000
STO-3G
```

---

## ⚙️ Running the Code

You can compile the code using a `Makefile`, then run using the provided `run.sh`:

```bash
bash run.sh
```

`run.sh` defines the input molecule, method, basis set, SCF settings, etc.

---

## 🚀 Features

- ✅ Supports **Restricted (RHF)** and **Unrestricted (UHF)** HF
- ✅ Configurable SCF options: **DIIS** or **Fixed-Point Iteration**
- ✅ Mulliken & Löwdin population analysis
- ✅ Modular: add new basis sets in `AOfunctions` in `mol.cpp`
- ✅ Output includes:
  - MO energies and coefficients
  - Convergence info and error tracking
  - Population analysis
  - Total energy

---

## 🧠 Libraries

| Library   | Description                                           |
|-----------|-------------------------------------------------------|
| `libmol`  | Parses input files and sets up atomic basis functions |
| `libmath` | Matrix operations, inversion, square roots, etc.     |
| `libgf`   | Gaussian primitives and Obara–Saika recursion         |
| `libint`  | One- and two-electron integral computation            |
| `libscf`  | SCF routines, DIIS, FPI for RHF/UHF                   |

---

## 📌 Sample Output

```
Total E     = -75.983948 Ha
Löwdin Pop. Analysis
idx              charge
1               -0.8343
2                0.4171
3                0.4171
Sum of atomic charges = 0.0000
```

---

## 🧪 Tested Molecules

| Molecule | Multiplicity | Method | HF Energy (Ha) |
|----------|---------------|--------|----------------|
| H₂       | 1             | RHF    | -1.11666       |
| H        | 2             | UHF    | -0.46658       |
| CH₃      | 2             | UHF    | -39.07670      |
| CH₄      | 1             | RHF    | -39.72674      |
| NH₃      | 1             | RHF    | -54.04445      |
| H₂O      | 1             | RHF    | -74.98394      |
| O₂       | 3             | UHF    | -144.09105     |

---

## 📜 License

MIT License © Prateek Saxena, BITS Pilani

> Because every electron deserves its own orbital...

# 📚 libscf: Self-Consistent Field (SCF) Module

This module handles all the SCF logic in your Hartree-Fock solver.

---

## 🔧 Files in `libscf`
- `scfgen.cpp` – Utility generators and helpers for SCF matrices and energy
- `scfrhf.cpp` – Restricted Hartree-Fock implementation
- `scfuhf.cpp` – Unrestricted Hartree-Fock implementation

---

## 📂 `scfgen.cpp`

### 🔹 `overlap(std::vector<GF> phis)`
**Return Type**: `Matrix`

**Description**: _Functionality to be filled based on function logic._

```cpp
overlap(std::vector<GF> phis){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			M.matrix[i][j] = S(phis[i],phis[j]);
		}
```

### 🔹 `kinetic(std::vector<GF> phis)`
**Return Type**: `Matrix`

**Description**: _Functionality to be filled based on function logic._

```cpp
kinetic(std::vector<GF> phis){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			M.matrix[i][j] = T(phis[i],phis[j]);
		}
```

### 🔹 `nuclear(std::vector<GF> phis, std::vector<int> Zvals, std::vector<std::vector<double>> xyzN)`
**Return Type**: `Matrix`

**Description**: _Functionality to be filled based on function logic._

```cpp
nuclear(std::vector<GF> phis, std::vector<int> Zvals, std::vector<std::vector<double>> xyzN){
	Matrix M(phis.size(),phis.size());
	for(int i = 0; i < phis.size(); i++){
		for(int j = 0; j < phis.size(); j++){
			for(int k = 0; k < Zvals.size(); k++){
				M.matrix[i][j] += -Zvals[k]*V(phis[i],phis[j],xyzN[k]);
			}
```

### 🔹 `nucrepl(std::vector<int> Z, std::vector<std::vector<double>> xyzN)`
**Return Type**: `double`

**Description**: _Functionality to be filled based on function logic._

```cpp
nucrepl(std::vector<int> Z, std::vector<std::vector<double>> xyzN){
	double sum = 0;
	double Rij;
	std::vector<double> Ri;
	std::vector<double> Rj;
	for(int i = 0; i < Z.size(); i++){
		for(int j = (i+1); j < Z.size(); j++){
			Ri = xyzN[i];
			Rj = xyzN[j];
			Rij = sqrt((Ri[0]-Rj[0])*(Ri[0]-Rj[0]) + 
				   (Ri[1]-Rj[1])*(Ri[1]-Rj[1]) + 
				   (Ri[2]-Rj[2])*(Ri[2]-Rj[2]));
			sum += Z[i]*Z[j] / Rij;
		}
```

## 📂 `scfrhf.cpp`

### 🔹 `R_density_matrix(Matrix C, int N)`
**Return Type**: `Matrix`

**Description**: _Functionality to be filled based on function logic._

```cpp
R_density_matrix(Matrix C, int N){
	Matrix P(C.rows, C.cols);
	for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			double sum = 0;
			for(int a = 0; a < N/2; a++){
				sum += C.matrix[i][a]*C.matrix[j][a];
			}
```

### 🔹 `R_F(Matrix Hcore, Matrix P, std::vector<std::vector<std::vector<std::vector<double>>>> g)`
**Return Type**: `Matrix`

**Description**: _Functionality to be filled based on function logic._

```cpp
R_F(Matrix Hcore, Matrix P, std::vector<std::vector<std::vector<std::vector<double>>>> g){
	Matrix G(P.rows, P.cols);
	double sum;
	for(int mu = 0; mu < G.rows; mu++){
		for(int nu = 0; nu < G.rows; nu++){
			sum = 0;
			for(int ld = 0; ld < G.rows; ld++){
				for(int sg = 0; sg < G.rows; sg++){
					sum += P.matrix[ld][sg] * (g[mu][nu][sg][ld] - 0.5*g[mu][ld][sg][nu]);
				}
```

### 🔹 `R_E0(Matrix P, Matrix Hcore, Matrix F)`
**Return Type**: `double`

**Description**: _Functionality to be filled based on function logic._

```cpp
R_E0(Matrix P, Matrix Hcore, Matrix F){
	double sum = 0;
	for(int i = 0; i < P.rows; i++){
		for(int j = 0; j < P.cols; j++){
			sum += P.matrix[j][i]*(Hcore.matrix[i][j]+F.matrix[i][j]);
		}
```

### 🔹 `R_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N)`
**Return Type**: `void`

**Description**: _Functionality to be filled based on function logic._

```cpp
R_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N){
	// Uses Ediff between iterations for error metric
	std::vector<Matrix> tec(2);
	double tEo;

	*f = R_F(hcore, *p, eris);
	tEo  = R_E0(*p, hcore, *f);
	*err = tEo - *Eo;
	*Eo = tEo;

	*fo  = transpose(x) * (*f) * x;
	tec  = diagonalize(*fo);
	*e   = tec[0];
	*co  = tec[1];
	*c   = x * (*co);
	*p   = R_density_matrix(*c, N);
}
```

### 🔹 `R_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps, int* icd)`
**Return Type**: `void`

**Description**: _Functionality to be filled based on function logic._

```cpp
R_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* p, Matrix* f, Matrix* fo, Matrix* e, Matrix* co, Matrix* c, double* Eo, double* err, int N, int i, Matrix* SPf, Matrix* SPe, int sps, int* icd){
	// Uses commutation of F and P for error metric
	// Perform fixed-point iterations until iteration number
	// equals the subspace size, then perform DIIS iterations
	if(i < sps){
		Matrix d = *p;
		R_FPI(s, hcore, eris, x, p, f, fo, e, co, c, Eo, err, N);
		for(int j = 0; j < sps-1; j++){
			SPf[j] = SPf[j+1];
			SPe[j] = SPe[j+1];
		}
```

### 🔹 `if(i < sps)`
**Return Type**: `then perform DIIS iterations`

**Description**: _Functionality to be filled based on function logic._

```cpp
if(i < sps){
		Matrix d = *p;
		R_FPI(s, hcore, eris, x, p, f, fo, e, co, c, Eo, err, N);
		for(int j = 0; j < sps-1; j++){
			SPf[j] = SPf[j+1];
			SPe[j] = SPe[j+1];
		}
```

### 🔹 `for(int j = 0; j < sps-1; j++)`
**Return Type**: `update error vectors`

**Description**: _Functionality to be filled based on function logic._

```cpp
for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			double sum = 0;
			for(int a = 0; a < N/2; a++){
				sum += C.matrix[i][a]*C.matrix[j][a];
			}
```

## 📂 `scfuhf.cpp`

### 🔹 `UR_density_matrix(Matrix C, int N)`
**Return Type**: `Matrix`

**Description**: _Functionality to be filled based on function logic._

```cpp
UR_density_matrix(Matrix C, int N){
	Matrix P(C.rows, C.cols);
	for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			double sum = 0;
			for(int a = 0; a < N; a++){
				sum += C.matrix[i][a]*C.matrix[j][a];
			}
```

### 🔹 `UR_F(Matrix Hcore, Matrix PT, Matrix Ps, std::vector<std::vector<std::vector<std::vector<double>>>> g)`
**Return Type**: `Matrix`

**Description**: _Functionality to be filled based on function logic._

```cpp
UR_F(Matrix Hcore, Matrix PT, Matrix Ps, std::vector<std::vector<std::vector<std::vector<double>>>> g){
	Matrix F(PT.rows, PT.cols);
	double sum;
	for(int mu = 0; mu < F.rows; mu++){
		for(int nu = 0; nu < F.rows; nu++){
			sum = Hcore.matrix[mu][nu];
			for(int ld = 0; ld < F.rows; ld++){
				for(int sg = 0; sg < F.rows; sg++){
					sum += (PT.matrix[ld][sg] * g[mu][nu][sg][ld]) - (Ps.matrix[ld][sg] * g[mu][ld][sg][nu]);
				}
```

### 🔹 `UR_E0(Matrix PT, Matrix Pa, Matrix Pb, Matrix Hcore, Matrix Fa, Matrix Fb)`
**Return Type**: `double`

**Description**: _Functionality to be filled based on function logic._

```cpp
UR_E0(Matrix PT, Matrix Pa, Matrix Pb, Matrix Hcore, Matrix Fa, Matrix Fb){
	double sum = 0;
	for(int i = 0; i < PT.rows; i++){
		for(int j = 0; j < PT.cols; j++){
			sum += PT.matrix[j][i] * Hcore.matrix[i][j] + Pa.matrix[j][i] * Fa.matrix[i][j] + Pb.matrix[j][i] * Fb.matrix[i][j];
		}
```

### 🔹 `UR_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb)`
**Return Type**: `void`

**Description**: _Functionality to be filled based on function logic._

```cpp
UR_FPI(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb){
	std::vector<Matrix> tec_a(2);
	std::vector<Matrix> tec_b(2);
	double tEo;

	*fa = UR_F(hcore, *pt, *pa, eris);
	*fb = UR_F(hcore, *pt, *pb, eris);
	
	tEo = UR_E0(*pt, *pa, *pb, hcore, *fa, *fb);
	*err = tEo - *Eo;
	*Eo = tEo;

	*fao  = transpose(x) * (*fa) * x;
	*fbo  = transpose(x) * (*fb) * x;
	tec_a  = diagonalize(*fao);
	tec_b  = diagonalize(*fbo);
	*ea  = tec_a[0];
	*cao = tec_a[1];
	*eb  = tec_b[0];
	*cbo = tec_b[1];
	*ca  = x * (*cao);
	*cb  = x * (*cbo);

	*pa  = UR_density_matrix(*ca, Na);
	*pb  = UR_density_matrix(*cb, Nb);
	*pt  = (*pa) + (*pb);
}
```

### 🔹 `UR_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, Matrix* SPfa, Matrix* SPfb, Matrix* SPea, Matrix* SPeb, int sps, int* icd)`
**Return Type**: `void`

**Description**: _Functionality to be filled based on function logic._

```cpp
UR_DIIS(Matrix s, Matrix hcore, std::vector<std::vector<std::vector<std::vector<double>>>> eris, Matrix x, Matrix* pt, Matrix* pa, Matrix* pb, Matrix* fa, Matrix* fb, Matrix* fao, Matrix* fbo, Matrix* ea, Matrix* eb, Matrix* cao, Matrix* cbo, Matrix* ca, Matrix* cb, double* Eo, double* err, int Na, int Nb, int i, Matrix* SPfa, Matrix* SPfb, Matrix* SPea, Matrix* SPeb, int sps, int* icd){
	// Uses commutation of F and P for error metric
	// Perform fixed-point iterations until iteration number
	// equals the subspace size, then perform DIIS iterations
	if(i < sps){
		Matrix da = *pa;
		Matrix db = *pb;
		UR_FPI(s, hcore, eris, x, pt, pa, pb, fa, fb, fao, fbo, ea, eb, cao, cbo, ca, cb, Eo, err, Na, Nb);
		for(int j = 0; j < sps-1; j++){
			SPfa[j] = SPfa[j+1];
			SPfb[j] = SPfb[j+1];
			SPea[j] = SPea[j+1];
			SPeb[j] = SPeb[j+1];
		}
```

### 🔹 `if(i < sps)`
**Return Type**: `then perform DIIS iterations`

**Description**: _Functionality to be filled based on function logic._

```cpp
if(i < sps){
		Matrix da = *pa;
		Matrix db = *pb;
		UR_FPI(s, hcore, eris, x, pt, pa, pb, fa, fb, fao, fbo, ea, eb, cao, cbo, ca, cb, Eo, err, Na, Nb);
		for(int j = 0; j < sps-1; j++){
			SPfa[j] = SPfa[j+1];
			SPfb[j] = SPfb[j+1];
			SPea[j] = SPea[j+1];
			SPeb[j] = SPeb[j+1];
		}
```

### 🔹 `for(int j = 0; j < sps-1; j++)`
**Return Type**: `update error vectors`

**Description**: _Functionality to be filled based on function logic._

```cpp
for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			double sum = 0;
			for(int a = 0; a < N; a++){
				sum += C.matrix[i][a]*C.matrix[j][a];
			}
```


# Numerical Methods in Python

![Python](https://img.shields.io/badge/language-Python-3776AB) ![NumPy](https://img.shields.io/badge/NumPy-SciPy-013243) ![Matplotlib](https://img.shields.io/badge/plots-Matplotlib-11557c)

*[Version française](README.fr.md)*

Six team projects from the **Numerical Algorithms** course (IS104) at [ENSEIRB-MATMECA](https://enseirb-matmeca.bordeaux-inp.fr/), 1st year, 2024–2025.
Each project implements numerical methods from scratch in Python, then applies them to a concrete problem.

| # | Project | Methods | Application |
|---|---------|---------|-------------|
| 1 | [Precision & CORDIC](01-precision-cordic) | Finite-precision arithmetic, CORDIC algorithms | Computing `ln`, `exp`, `arctan`, `tan` |
| 2 | [Linear systems](02-cholesky-equation-chaleur) | Full / incomplete Cholesky, preconditioned conjugate gradient | Heat equation |
| 3 | [SVD image compression](03-svd-compression-image) | Householder, bidiagonalisation, QR, SVD | Lossy image compression |
| 4 | [Non-linear systems](04-newton-raphson) | Newton–Raphson with backtracking | Lagrange points, electrostatic equilibrium |
| 5 | [Interpolation & integration](05-splines-integration) | Cubic splines, rectangle, trapezoid, Simpson rules | Airfoil profile and lift |
| 6 | [Differential equations](06-equations-differentielles) | Euler, midpoint, Heun, Runge–Kutta 4 | Double pendulum (chaos), predator–prey model |

Projects 2, 4, 5 and 6 include a report (`rapport.pdf`, in French) describing the approach and results.

<p align="center">
  <img src="05-splines-integration/rapport/airflow.png" width="45%" alt="Airflow around the airfoil">
  <img src="06-equations-differentielles/rapport/img/extremity_double_pendulum.png" width="45%" alt="Trajectory of the tip of a double pendulum">
</p>

## Installation

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Projects

### 1. Precision and CORDIC algorithms
Simulates a machine working with `p` significant digits (`reduce`) and studies relative errors on addition and multiplication. Computes `ln`, `exp`, `arctan` and `tan` with CORDIC-style methods (as in a pocket calculator). `projet1bis.py` is the final version, with the error analysis.

```bash
cd 01-precision-cordic && python3 projet1bis.py
```

### 2. Solving linear systems — heat equation
Full and incomplete Cholesky factorisation, conjugate gradient with and without preconditioning, applied to the 2D heat equation.

```bash
cd 02-cholesky-equation-chaleur && python3 test.py
```

### 3. SVD image compression
Householder transformations, reduction to bidiagonal form, QR iterations to compute the SVD, then compression of an RGB image by truncating singular values.

```bash
cd 03-svd-compression-image/src && python3 bidiagonale.py
cd ../partie-4 && python3 traitement_image.py
```

### 4. Non-linear systems — Newton–Raphson
Newton–Raphson with backtracking, applied to the 5 Lagrange points of a two-body system and to the electrostatic equilibrium of `N` charges (related to the roots of Legendre polynomials).

```bash
cd 04-newton-raphson && python3 lagrange.py && python3 EquilibreES.py
```

### 5. Interpolation and integration — airfoil
Interpolates an airfoil profile (`e662.dat`) with cubic splines, compares the convergence speed of integration methods, then computes the pressure map and the lift.

```bash
cd 05-splines-integration/src && python3 cl_airfoil_analysis.py e662.dat
```

### 6. Differential equations — dynamical systems
One-step methods (Euler, midpoint, Heun, RK4) with fixed or adaptive step size. Studies the simple pendulum, the chaotic double pendulum, and the Malthus, Verhulst and Lotka–Volterra models.

```bash
cd 06-equations-differentielles && python3 src/pendulum.py && python3 src/proie_predateur.py
```

## Authors

Each project was done in a team of 3 or 4 students:

| Project | Team |
|---------|------|
| 1 | Yani Haouchine, Antoine Glasson, Benjamin Chew, Ihab Jaouad |
| 2 | Yani Haouchine, Arnaud Chane-Kive, Wissal Faouzi, Vivien Nebout |
| 3 | Yani Haouchine, Benjamin Chew, Tim Lacault |
| 4 | Yani Haouchine, Benjamin Basset, Margaux Lopez, Yoann Moresco |
| 5 | Yani Haouchine, Yoann Moresco, Wissal Faouzi, Solan Boher-Sanchez |
| 6 | Yani Haouchine, Benjamin Chew, Solan Boher-Sanchez, Alexis Gaudray Bouju |

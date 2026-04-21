# Naive Gauss-Jordan Solver

A lightweight, pure-Python implementation of the Gauss-Jordan elimination method for linear algebra operations. This library is designed for educational purposes and for scenarios where **exact rational arithmetic** is required.

## ✨ Features

- **Exact Arithmetic Support**: By using Python's `fractions.Fraction` objects, you can perform calculations without floating-point rounding errors, obtaining mathematically perfect results.
- **Linear Equation Solver**: Solve systems of linear equations ($Ax = b$).
- **Matrix Inversion**: Compute the inverse of a square matrix.
- **Matrix Operations**: 
  - Matrix Multiplication (Dot product)
  - Matrix Transposition
  - Determinant calculation
  - Matrix Rank calculation
  - Moore-Penrose Pseudoinverse
- **Zero Dependencies**: Built using only the Python Standard Library. No `numpy` required.
- **Pretty Printing**: Built-in utilities to display matrices in a clean, readable format.

## 🚀 Installation

Since this is a single-file library, you can simply download `gaussjordan.py` and import it into your project.

```bash
# Download the file
curl -O https://raw.githubusercontent.com/kitanokitsune/naive_gauss_jordan/main/gaussjordan.py
```

## 🛠 Usage

### 1. Solving Linear Equations
Solve the system:  
$1x + 2y + 5z = 12$  
$1x + 3y + 1z = 10$  
$2x + 3y + 1z = 14$  

```python
from gaussjordan import solve

M = [[1, 2, 5], 
     [1, 3, 1], 
     [2, 3, 1]]
V = [12, 10, 14]

ans = solve(M, V)
print(f"Solution: {ans}") 
# Output: [4.0, 1.6923076923076916, 0.9230769230769231]
```

### 2. Exact Arithmetic with Fractions
To avoid floating-point errors, convert your matrices to `Fraction` objects before processing.

```python
from gaussjordan import solve, toFraction

M = [[1, 2, 5], [1, 3, 1], [2, 3, 1]]
V = [12, 10, 14]

# Convert to exact rational numbers
Mf = toFraction(M)
Vf = toFraction(V)

# Solve exactly
ans_exact = solve(Mf, Vf)

print(f"Exact Solution: {ans_exact}")
# Output: [Fraction(4, 1), Fraction(22, 13), Fraction(12, 13)]
```

### 3. Matrix Inversion and Pretty Printing
```python
from gaussjordan import invert, toFraction, pprint_mat

M = toFraction([ 
     [1, 2], 
     [3, 4]])

inv_M = invert(M)
pprint_mat("Inverse Matrix=", inv_M)
# Output: 
# Inverse Matrix={  -2,    1 |
#                | 3/2, -1/2 }
```

### 4. Moore-Penrose Pseudoinverse
```python
from gaussjordan import rank, moore_penrose, toFraction, pprint_mat

M = toFraction([
     [1, 2, 3, 1],
     [2, 4, 6, 4],
     [0, 0, 0, 2]])

r = rank(M)
print(f"Rank: {r}")
# Output: 2

pinv_M = moore_penrose(M)
print("Moore-Penrose Pseudoinverse:")
pprint_mat("", pinv_M)
# Output: 
# { 1/28, 1/56, -3/56 |
# | 1/14, 1/28, -3/28 |
# | 3/28, 3/56, -9/56 |
# | -1/6, 1/12,  5/12 }
```


## 📖 API Reference

| Function | Description |
| :--- | :--- |
| `solve(A, b)` | Solves the linear system $Ax = b$. |
| `invert(mat)` | Returns the inverse of a square matrix. |
| `dot(A, B)` | Performs matrix multiplication (dot product). |
| `det(A)` | Calculates the determinant of a matrix. |
| `rank(mat)` | Returns the rank of the matrix. |
| `moore_penrose(mat)` | Calculates the Moore-Penrose pseudoinverse. |
| `transpose(A)` | Returns the transpose of a matrix. |
| `toFraction(M)` | Recursively converts a matrix to `fractions.Fraction`. |
| `toReal(M)` | Recursively converts elements to `float`. |
| `pprint_mat(prefix="", mat=[])` | Prints the matrix in a formatted, readable way. |
| `spprint_mat(prefix="", mat=[])` | Returns the output of `pprint_mat()` as a string. |

## ⚠️ Disclaimer

This implementation is a **"Naive"** solver. While it is excellent for learning and exact rational arithmetic, it is not optimized for extremely large matrices or high-performance computing like `NumPy`.

## 📜 License

This project is released under the [MIT License](LICENSE).


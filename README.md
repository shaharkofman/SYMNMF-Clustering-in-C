## 💻 README: Symmetric Non-negative Matrix Factorization (SymNMF) for Clustering

This repository contains the implementation of a clustering algorithm based on **Symmetric Non-negative Matrix Factorization (SymNMF)**, along with a comparative analysis against the **K-means** algorithm.

The core logic is implemented in C for performance, with a Python interface for user interaction and analysis.

-----

### 📝 Project Overview

The project implements the SymNMF clustering algorithm, which aims to find a non-negative decomposition matrix $\mathbf{H}$ by minimizing $min_{H\ge0}||W-HH^{T}||_{F}^{2}$, where $W$ is the Normalized Similarity Matrix (Graph Laplacian).

The full SymNMF process includes:

1.  Forming the **Similarity Matrix ($A$)**.
2.  Computing the **Diagonal Degree Matrix ($D$)**.
3.  Computing the **Normalized Similarity Matrix ($W$)**.
4.  Finding $\mathbf{H}$ using an iterative update rule until convergence or max iterations are reached.
5.  Deriving a **hard clustering solution** by assigning each data point to the cluster (column of $\mathbf{H}$) with the highest association score.

-----

### 📦 File Structure

The project consists of the following required files:

| File Name | Description |
| :--- | :--- |
| **`symnmf.py`** | Python interface for reading arguments, handling **H initialization**, and calling the C extension functions (`symnmf`, `sym`, `ddg`, `norm`). |
| **`symnmf.c`** | C implementation of the core mathematical functions and the full SymNMF iteration logic. Also supports command-line execution for `sym`, `ddg`, and `norm` goals. |
| **`symnmf.h`** | C header file defining function prototypes used by `symnmf.c` and `symnmfmodule.c`. |
| **`symnmfmodule.c`** | Python C API wrapper defining the C extension functions for use in Python. |
| **`analysis.py`** | Program to compare SymNMF clustering against **K-means** and report the **`silhouette_score`**. |
| **`setup.py`** | Build script used by Python to create the shared object (`.so`) file for the C extension. |
| **`Makefile`** | Script to build the standalone C executable (`./symnmf`). |

-----

### ⚙️ Build and Run Instructions

The project can be run via the Python interface (using the C extension) or as a standalone C executable.

#### 1\. Building the Python C Extension

This command builds the C extension, making the `symnmf` module available for import in Python.

```bash
python3 setup.py build_ext --inplace
```

#### 2\. Building the Standalone C Executable

The `Makefile` compiles the C program into an executable named `./symnmf`.

```bash
make
```

-----

### 🚀 Execution

#### 1\. Python Interface (`symnmf.py`)

Handles four goals: `symnmf` (full clustering), `sym`, `ddg`, or `norm`.

**Usage:**

```bash
python3 symnmf.py <k> <goal> <file_name.txt>
```

| Goal | Output |
| :--- | :--- |
| **`symnmf`** | The final decomposition matrix **H**. |
| **`sym`** | The **Similarity Matrix (A)**. |
| **`ddg`** | The **Diagonal Degree Matrix (D)**. |
| **`norm`** | The **Normalized Similarity Matrix (W)**. |

**Example:**

```bash
python3 symnmf.py 2 symnmf input_1.txt
```

#### 2\. C Standalone Program (`./symnmf`)

Supports three goals: `sym`, `ddg`, or `norm`.

**Usage:**

```bash
./symnmf <goal> <file_name.txt>
```

**Example:**

```bash
./symnmf sym input_1.txt
```

#### 3\. Analysis (`analysis.py`)

Compares SymNMF and K-means clustering performance using the silhouette score.

**Usage:**

```bash
python3 analysis.py <k> <file_name.txt>
```

**Example Output Format:**

```
nmf: 0.1162
kmeans: 0.1147
```

-----

### ⚠️ Implementation Details

  * **Output Format:** All matrix outputs must be formatted to **4 decimal places** (`%.4f`), with values separated by commas, and each row on a new line.
  * **H Initialization:** In `symnmf.py`, initialize H using `np.random.seed(1234)` and `np.random.uniform()` with values from $[0, 2\sqrt{m/k}]$, where $m$ is the average of all entries in $W$.
  * **Convergence Parameters:** Use $max\_iter=300$ and $\epsilon=1e-4$ for convergence of both K-means and NMF.
  * **Error Handling:** In case of any error, print `"An Error Has Occurred"` and terminate. All dynamically allocated memory must be freed.

This repository contains the implementation of a clustering algorithm based on **Symmetric Non-negative Matrix Factorization (SymNMF)**, along with a comparative analysis against the **K-means** algorithm.

The core logic is implemented in C for performance, with a Python interface for user interaction and analysis.

### 📝 Project Overview

The project implements the SymNMF clustering algorithm, which involves the following steps:

1.  [cite\_start]Calculate the **Similarity Matrix ($A$)**[cite: 9, 11].
2.  [cite\_start]Calculate the **Diagonal Degree Matrix ($D$)**[cite: 9, 16].
3.  [cite\_start]Calculate the **Normalized Similarity Matrix ($W$)** (Graph Laplacian)[cite: 9, 23, 25].
4.  [cite\_start]Find the non-negative decomposition matrix $\mathbf{H}$ by minimizing $min_{H\ge0}||W-HH^{T}||_{F}^{2}$ using an iterative update rule[cite: 9, 27, 28, 33].
5.  [cite\_start]Derive the final **hard clustering solution** from $\mathbf{H}$ by selecting the cluster with the highest association score for each data point[cite: 37, 39].

[cite\_start]The implementation is split between a C extension and Python files to manage different aspects of the assignment, including a required analysis against K-means using the Silhouette Score[cite: 112].

-----

### 📦 File Structure

The project consists of the following required files:

| File Name | Description |
| :--- | :--- |
| **`symnmf.py`** | [cite\_start]Python interface for reading arguments, handling **SymNMF initialization** [cite: 68][cite\_start], and calling the C extension functions (`symnmf`, `sym`, `ddg`, `norm`)[cite: 59, 73, 74, 75, 77]. |
| **`symnmf.c`** | [cite\_start]C implementation of the core mathematical functions and the full SymNMF iteration logic[cite: 53, 86]. [cite\_start]Also supports command-line execution for `sym`, `ddg`, and `norm` goals[cite: 89, 90, 91, 92]. |
| **`symnmf.h`** | [cite\_start]C header file defining function prototypes used by `symnmf.c` and `symnmfmodule.c`[cite: 109, 110]. |
| **`symnmfmodule.c`** | [cite\_start]Python C API wrapper defining the C extension functions (`symnmf`, `sym`, `ddg`, `norm`) for use in Python[cite: 104, 108]. |
| **`analysis.py`** | [cite\_start]Program to compare SymNMF clustering against **K-means** using the **`silhouette_score`** from `sklearn.metrics`[cite: 112, 113]. |
| **`setup.py`** | [cite\_start]Build script used by Python to create the shared object (`.so`) file for the C extension[cite: 56, 129]. |
| **`Makefile`** | [cite\_start]Script to build the standalone C executable (`./symnmf`)[cite: 57, 131]. |

-----

### ⚙️ Build and Run Instructions

The project has two main ways to run: via the Python interface (using the C extension) and as a standalone C executable for the matrix calculations.

#### 1\. Building the Python C Extension

This command builds the C extension, making the `symnmf` module available for import in Python (`symnmf.py` and `analysis.py`).

```bash
python3 setup.py build_ext --inplace
```

#### 2\. Building the Standalone C Executable

The `Makefile` is used to compile the C program (`symnmf.c` and its dependencies) into an executable named `./symnmf`. [cite\_start]Compilation uses the flags: `-ansi -Wall -Wextra -Werror -pedantic-errors`[cite: 132, 133].

```bash
make
```

-----

### 🚀 Execution

#### 1\. Python Interface (`symnmf.py`)

[cite\_start]The Python program handles four possible goals: `symnmf`, `sym`, `ddg`, or `norm`[cite: 62, 63, 64, 65, 66].

**Usage:**

```bash
python3 symnmf.py <k> <goal> <file_name>
```

**Examples:**

  * **Full SymNMF:** Outputs the final matrix $\mathbf{H}$.
    ```bash
    python3 symnmf.py 2 symnmf input_1.txt
    ```
  * **Similarity Matrix:** Outputs the similarity matrix $\mathbf{A}$.
    ```bash
    python3 symnmf.py 2 sym input_1.txt
    ```

#### 2\. C Standalone Program (`./symnmf`)

[cite\_start]The C program supports three goals: `sym`, `ddg`, or `norm`[cite: 89, 90, 91, 92].

**Usage:**

```bash
./symnmf <goal> <file_name>
```

**Example:**

  * **Similarity Matrix:**
    ```bash
    ./symnmf sym input_1.txt
    ```

#### 3\. Analysis (`analysis.py`)

[cite\_start]This program compares the SymNMF clustering result against K-means and reports the `silhouette_score` for both[cite: 112].

**Usage:**

```bash
python3 analysis.py <k> <file_name>
```

**Example:**

```bash
python3 analysis.py 5 input_k5_d7.txt
# Expected Output (example format):
# nmf: 0.1162
# kmeans: 0.1147
```

-----

### ⚠️ Important Notes and Assumptions

  * [cite\_start]**Floating Point Format:** All matrix outputs (in both C and Python) must be formatted to **4 decimal places** (`%.4f`)[cite: 143].
  * [cite\_start]**Data Types:** Use `double` in C and `float` in Python for vector/matrix elements[cite: 153].
  * [cite\_start]**Convergence:** For both K-means and NMF, use $max\_iter=300$ and $\epsilon=1e-4$[cite: 154].
  * [cite\_start]**Error Handling:** In case of any error, the program must print `"An Error Has Occurred"` and terminate[cite: 150]. [cite\_start]Don't forget to **free all allocated memory**[cite: 151].
  * [cite\_start]**No Argument Validation:** It is assumed that command-line arguments are legal[cite: 142].

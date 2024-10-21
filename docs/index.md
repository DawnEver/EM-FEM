
# Docs
## Geometry and Mesh
### Read Mesh from Gmsh
Tutorial for Gmsh can be found in [here](http://www.gmsh.info/doc/texinfo/gmsh.html).

When modelling with Gmsh, boundary conditions can be defined with Physical Line, and material can be define with Physical Surface,
- `Physical` section in `.geo file`(e.g. [test-T_A_Lipo.geo](./data//lipo/test-T_A_Lipo.geo)) to define information in the name of `Physical Surface` or `Physical Line`.
- Function `read_gmsh()` in [`read_mesh.py`](./src/read_mesh.py) to get definitions from string.

### Read Mesh from other files
When reading mesh from .nas(Nastran) or .csv with only vertex and triangle information, additional definitions of boundary conditions, group and material are required.
- Refer to function `read_nastran()` and `read_lipo_csv()` in [`read_mesh.py`](./src/read_mesh.py).

## Current Excitation
The current excitation is determined based on the group ID.
Rather than predefining a constant copper area, the total area of all triangles within a group is summed up to calculate the copper area.
Then is the calculation of current density for each group.

## Material Database
Refer to [`material.py`](./src/material.py).

Dict `material_db_dict` define all materials available here. The database can be read from an external file or database.

In one FEA case, program don't need to load all materials. To address this, Class `Material_In_Use` has been implemented to instantiate each material using the `Material` class. The `Material` class contains standardized function `reluctivity_func` to calculate reluctivity, with reluctivity for steel determined by BH curve interpolation and for other materials by relative permeability.

## Magneto-static Solver
The problem can finally be abstracted as solving linear systems:

$$
[S]⋅ [A] = [T]
$$

The reluctivity matrix $[S]$ is calculated with reluctivity and the geometric characteristics of the vertices of each triangle.
The current matrix $[T]$ is obtained by simply taking 1/3 the sum of the currents passing through the triangles.
More details can be found in *chapter 10, [Introduction to AC machine design]((https://onlinelibrary.wiley.com/doi/book/10.1002/9781119352181))*.

The magnetic vector potential matrix $[A]$ can be solved by direct method and iteration methods.

The NumPy package in Python offers a direct solving method for linear systems. To understand how it works, we can trace down its implementation. NumPy utilizes [LAPACK (Linear Algebra Package)](https://netlib.org/lapack/), which is written in Fortran, to solve linear systems using LU decomposition(Cholesky decomposition for positive definite matrix), forward substitution, and backward substitution.

Similarly, SciPy sparse package uses UMFPACK(included in [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse), a suite of sparse matrix algorithms) to solve unsymmetric sparse linear systems.

Besides, nonlinear material can provide different reluctivity at each epoch. After a certain number of iterations, $[A]$ will eventually converge below the required error.

## Solve large-scale linear systems
Sometimes, the linear system for extra-fine mesh is too extensive to solve easily and quickly.
There are several possible solutions.
### 1. Apply iteration methods like gauss-seidel method to get an approximate solution.
But with limited time and space the problem still can not converge very well. Advanced methods like successive over relaxation(SOR) iteration and conjugate gradient method also may not significantly enhance the convergence rate.

### 2. Assign proper initial value with interpolation functions to accelerate convergence.
While this approach holds theoretical promise, it raises a new practical challenge of determining the suitable initial value.

### 3. Use sparse matrix to reduce space and time requirements
In systems of finite element problems , most of the coefficients are typically zero.
Using sparse matrix allows for more efficient storage, faster computations, reduced memory requirements.

### 4. Rewrite with high performance programming language(C++ or Rust).
According to Knuth's optimization principle, "Premature optimization is the root of all evil."
In the context of a finite element problem, there are two time-consuming steps: preparing matrices and solving the system.
Python, when solving linear systems with libraries like NumPy or SciPy, internally calls established math libraries.
While rewriting the program in C++ might speed up the process of preparing matrices, it would not significantly enhance the solving performance.

## Post Processing
Energy, magnetic flux density can be calculated after obtaining magnetic vector potential of each vertex.

## Plot
Refer to [`plot.py`](./src/plot.py).

The supported plot types include scatter, counter, and vector plots, as defined in the Enum `PlotMapType`.
Interpolation also used to generate smoother counter plots.


## Test Cases
- [Simple Iron-cored Inductor with an Air Gap](./test_inductor.md)
- [A Synchronous Reluctance Machine](./test_synrm.md)

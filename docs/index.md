
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
The problem can finally be abstracted as solving linear equations:

$$
[S]⋅ [A] = [T]
$$

The reluctivity matrix $[S]$ is calculated with reluctivity and the geometric characteristics of the vertices of each triangle.
The current matrix $[T]$ is obtained by simply taking 1/3 the sum of the currents passing through the triangles.
More details can be found in *chapter 10, [Introduction to AC machine design]((https://onlinelibrary.wiley.com/doi/book/10.1002/9781119352181))*.

The magnetic vector potential matrix $[A]$ can be solved by gaussian elimination(more precise but slower) and gauss-seidel iteration(large quantity of equations).

Besides, nonlinear material can provide different reluctivity at each epoch. After a certain number of iterations, $[A]$ will eventually converge below the required error.

## Post Processing
Energy, magnetic flux density can be calculated after obtaining magnetic vector potential of each vertex.

## Plot
Refer to [`plot.py`](./src/plot.py).

The supported plot types include scatter, counter, and vector plots, as defined in the Enum `PlotMapType`.
Interpolation also used to generate smoother counter plots.


## Test Cases
- [Simple Iron-cored Inductor with an Air Gap](./test_lipo.md)
- [A Synchronous Reluctance Machine](./test_synrm.md)

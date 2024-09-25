# Electromagnetic Finite Element Method
Please read [Docs](./docs/index.md) for more details.

## Dependencies
- [gmsh](https://gmsh.info)
- python>=3.10
- [requirements.txt](./requirements.txt)

## Quick Start
1. Install dependencies.
```sh
pip install -r requirements.txt
```
2. Generate mesh with Gmsh.
```sh
gmsh
# or generate from .geo file
gmsh <path/to/file.geo>
```
3. Define current excitation, solve and plot in python script. Below are 2 test cases:
- [test_lipo.py](./src/test_lipo.py) Test case of a simple iron-cored inductor with an air gap.
    - coarse mesh same as Lipo's example(.csv)
    - extra fine mesh generated with gmsh(.msh)
- [test_synrm.py](./src/test_synrm.py) Test case of a synchronous reluctance machine.
    - coarse mesh(.nas)
    - extra fine mesh(.nas)

## Docs
- [docs](./docs/index.md)
- [test case 1](./docs/test_lipo.md)
- [test case 2](./docs/test_synrm.md)

## Dev
1. Install dependencies.
```sh
pip install -r requirements-dev.txt
```

2. [Ruff](https://github.com/astral-sh/ruff) to format python scripts.

3. Install pre-commit git hooks in [.pre-commit-config.yaml](./.pre-commit-config.yaml)
```sh
pre-commit install .
```

4. [Makefile](./Makefile) defines some useful tools.

## TODO
- rewrite solver with high performance programming language(C++ or Rust).
- Assign proper initial value with interpolation functions to accelerate convergence.
- time step simulation & external circuit excitation
- PM materials
- parametric modelling of several classical electric machines
- post processing like loss calculation(Steinmetz model), torque(Maxwell stress tensor method)
    - loss calculation methods
        - calculate with maximum flux density(quick but coarse)
        - calculate with the amplitude of each harmonic flux density(more precise)

## Ref
- [T. A. Lipo, Introduction to AC machine design. in IEEE Press series on power engineering. Hoboken, New Jersey: IEEE Press/Wiley, 2017.](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119352181)
- [Gitbook Numerical Analysis Implementations in Python](https://laus-organization.gitbook.io/numerical-analysis-implementations-in-python/)

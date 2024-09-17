# Test Case: A Synchronous Reluctance Machine(SynRM)
Refer to [test_synrm.py](../src/test_synrm.py)

<img src='./test_synrm/synrm_model.png' width=50%>

[Test case of inductor](./test_lipo.md) has proven the accuracy and extensibility of my code, my code can be applied to the model of SynRM with coarse and extra fine mesh, while only modifying mesh and current excitation.

For large-scale systems of linear equations, it may be difficult for classical methods like gaussian elimination and LU decomposition to solve.
Iteration methods like gauss-seidel method are developed to get an approximate solution.

As table is shown below, for extra fine mesh of SynRM in this case, the system of equations of is too large to solve with classical method, gauss-seidel method is used.

| | Coarse Mesh | Extra Fine Mesh |
|--|--|--|
| Total Flux(Wb/m) | 0.712218 | --  |
| Total Energy(W) | 1121.6359 | --  |
| Num. of Vertex | 5540 | 85912  |

To be honest, the FEA of extra fine mesh still can not converge very well(only tried with limited iterations).
We can improve this in two ways:
1. Try other advanced methods like successive over relaxation(SOR) iteration and conjugate gradient method.
2. Rewrite the solver component with high performance programming language(C++ or Rust).

Work is still ongoing.

<table>
    <tr>
      <th><img src='./test_synrm/coarse-A.png' height=300pt><br/>Magnetic vector potential(Coarse)</th>
      <th><img src='./test_synrm/fine-A.png' height=300pt><br/>Magnetic vector potential(Fine)</th>
    </tr>
    <tr>
      <th><img src='./test_synrm/coarse-B.png' height=300pt><br/>FLux density(Coarse)</th>
      <th><img src='./test_synrm/fine-B.png' height=300pt><br/>FLux density(Fine)</th>
    </tr>
    <tr>
      <th><img src='./test_synrm/coarse-B_norm.png' height=300pt><br/>FLux density(Coarse)</th>
      <th><img src='./test_synrm/fine-B_norm.png' height=300pt><br/>FLux density(Fine)</th>
    </tr>
    <tr>
      <th><img src='./test_synrm/coarse-Energy.png' height=300pt><br/>Energy(Coarse)</th>
      <th><img src='./test_synrm/fine-Energy.png' height=300pt><br/>Energy(Fine)</th>
    </tr>
    <tr>
      <th><img src='./test_synrm/coarse-J.png' height=300pt><br/>Current density(Coarse)</th>
      <th><img src='./test_synrm/fine-J.png' height=300pt><br/>Current density(Fine)</th>
    </tr>
</table>

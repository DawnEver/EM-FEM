# Test Case: Simple Iron-cored Inductor with an Air Gap
Refer to [test_lipo.py](../src/test_lipo.py)


<table>
    <tr>
      <th><img src='./test_lipo/inductor_model_lipo.jpg' height=300pt><br/>Mesh in Lipo's book.</th>
      <th><img src='./test_lipo/inductor_model_gmsh.png' height=300pt><br/>Mesh in Gmsh</th>
    </tr>
</table>

Starting from a simple iron-cored inductor with an air gap discussed in *[Introduction to AC machine design]((https://onlinelibrary.wiley.com/doi/book/10.1002/9781119352181))*, T. A. Lipo provides example code in the appendix. (Refer to [Lipo.m](./test_lipo/Lipo.m))

The example code provide a simple practice but have several shortcomings:
- Low-quality mesh: For the simplicity of counter plot(Flux density plot), which needs grid-like points, the triangles are not optimal.
- Low scalability: Modeling, meshing, solving, and post-processing are highly coupled in the example code, which limits the program's extensibility.
- Too much repetitive code ...

In my code, mesh, Gmsh is used to model.
Mesh, solver and post-processing components are decoupled.
User can import and analyze any Gmsh geometry with only few changes(mesh and current excitation).

## Results Comparison
To verify the accuracy of the results of my python code.
I export the coarse mesh in Lipo's code and use Gmsh to generate finer mesh.
Then I compared the FEA results of Lipo's code and my code with the same mesh(Lipo's code not support other mesh).

| | Lipo's Code | Lipo's Mesh + My code | Gmsh Mesh + My code |
|--|--|--|--|
| Total Flux(Wb/m) | 0.0149 | 0.014908  | 0.014290 |
| Total Energy(J) | 28.1353 | 28.1353  | 27.0102 |
| Num. of Vertex | 1428 | 1428  | 6466 |

<table>
    <tr>
      <th><img src='./test_lipo/A-lipo.jpg' height=300pt><br/>Magnetic vector potential plot of Lipo's Code</th>
      <th><img src='./test_lipo/A-py-coutour-coarse.png' height=300pt><br/>Magnetic vector potential plot of my Code</th>
    </tr>
</table>

---
Below are some other plots obtained from my code:

<table>
    <tr>
      <th><img src='./test_lipo/A-py-coutour.png' height=300pt><br/>Magnetic vector potential</th>
      <th><img src='./test_lipo/A-py-coutourf.png' height=300pt><br/>Magnetic vector potential</th>
    </tr>
    <tr>
      <th><img src='./test_lipo/B-py.png' height=300pt><br/>FLux density</th>
      <th><img src='./test_lipo/B_norm-py.png' height=300pt><br/>FLux density</th>
    </tr>
    <tr>
      <th><img src='./test_lipo/J-py.png' height=300pt><br/>Current density</th>
      <th><img src='./test_lipo/Energy-py.png' height=300pt><br/>Energy</th>
    </tr>
</table>

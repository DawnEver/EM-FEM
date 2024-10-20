// test symmetric boundary condition

Geometry.AutoCoherence = 0;
// Define mesh

mesh_size_demo_steel = 3;
mesh_size_cu = 3;

// define points
Point(1) = {0, 0, 0, mesh_size_demo_steel};
Point(2) = {6, 0, 0, mesh_size_demo_steel};
Point(3) = {6, 4, 0, mesh_size_demo_steel};
Point(4) = {0, 4, 0, mesh_size_demo_steel};

Point(5) = {4, 1, 0, mesh_size_cu};
Point(6) = {4, 3, 0, mesh_size_cu};
Point(7) = {5.5, 3, 0, mesh_size_cu};
Point(8) = {5.5, 1, 0, mesh_size_cu};

// points to lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// generate surface
Coherence;
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1,2};
Plane Surface(2) = {2};


// physical surface
Physical Surface("group@s1@steel@demo_steel") = {1};
Color Blue {Surface{1};}
Physical Surface("group@c1@copper@cu") = {2};
Color Orange {Surface{2};}

// boundary condition examples
Physical Line("boundary@dirichlet@0") = {1,3};

Physical Line("boundary@symmetry@even1") = {4};
Physical Line("boundary@symmetry@even2") = {2};

Mesh 2;
Save "test_sym_bc.msh";

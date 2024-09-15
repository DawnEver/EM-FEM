Geometry.AutoCoherence = 0;
// Define mesh

mesh_size_demo_steel = 3;
mesh_size_cu = 3;
mesh_size_air = 3;

// define points
Point(1) = {-50, 50, 0, mesh_size_demo_steel};
Point(2) = {50,-50, 0, mesh_size_demo_steel};
Point(3) = {-50,-50, 0, mesh_size_demo_steel};
Point(4) = {50,50, 0, mesh_size_demo_steel};
Point(5) = {-10,10, 0, mesh_size_cu};
Point(6) = {10,-10, 0, mesh_size_cu};
Point(7) = {-10,-10, 0, mesh_size_cu};
Point(8) = {10,10, 0, mesh_size_cu};
Point(9) = {-30,-30, 0, mesh_size_air};
Point(10) = {-30,-40, 0, mesh_size_air};
Point(11) = {30,-40, 0, mesh_size_air};
Point(12) = {30,-30, 0, mesh_size_air};


// points to lines
Line(1) = {1, 3};
Line(2) = {3, 2};
Line(3) = {2, 4};
Line(4) = {4, 1};
Line(5) = {5, 7};
Line(6) = {7, 6};
Line(7) = {6, 8};
Line(8) = {8, 5};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

// generate surface
Coherence;
Curve Loop(1) = {8, 5, 6, 7};
Curve Loop(2) = {4, 1, 2, 3};
Curve Loop(3) = {12, 9, 10, 11};
Plane Surface(1) = {1};
Plane Surface(2) = {1,2,3};
Plane Surface(3) = {3};

// physical surface
Physical Surface("group@c1@copper@cu") = {1};
Color Red {Surface{1};}
Physical Surface("group@s1@steel@demo_steel") = {2};
Color Green {Surface{2};}
Physical Surface("group@a1@other@air") = {3};
Color Orange {Surface{3};}

// boundary condition examples
// Physical Line("boundary@dirichlet@0") = {4};
// Physical Line("boundary@dirichlet@1e-2") = {2};
// Physical Line("boundary@symmetry@odd1") = {2};
// Physical Line("boundary@symmetry@odd2") = {4};
// Physical Line("boundary@dirichlet@0") = {1,2,3,4};

Mesh 2;
Save "test-T_A_Lipo.msh";

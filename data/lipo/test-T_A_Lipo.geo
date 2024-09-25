Geometry.AutoCoherence = 0;
// Define mesh


mesh_size_demo_steel = 0.4;
mesh_size_air = 0.1;

// define points
Point(1) = {0, 0, 0, mesh_size_demo_steel};
Point(2) = {6, 0, 0, mesh_size_demo_steel};
Point(3) = {6, 2, 0, mesh_size_air};
Point(4) = {0, 2, 0, mesh_size_air};

Point(5) = {0, 2.5, 0, mesh_size_air};
Point(6) = {2, 2.5, 0, mesh_size_air};
Point(7) = {2, 4.5, 0, mesh_size_air};
Point(8) = {4, 4.5, 0, mesh_size_air};
Point(9) = {4, 2.5, 0, mesh_size_air};
Point(10) = {6, 2.5, 0, mesh_size_air};
Point(11) = {6, 6.5, 0, mesh_size_demo_steel};
Point(12) = {0, 6.5, 0, mesh_size_demo_steel};

Point(13) = {0-1, 0-0.1, 0, mesh_size_air};
Point(14) = {6+1, 0-0.1, 0, mesh_size_air};
Point(15) = {6+1, 6.5+0.1, 0, mesh_size_air};
Point(16) = {0-1, 6.5+0.1, 0, mesh_size_air};

Point(17) = {2+0.1, 2.5, 0, mesh_size_air};
Point(18) = {4-0.1, 2.5, 0, mesh_size_air};
Point(19) = {4-0.1, 4.5-0.1, 0, mesh_size_air};
Point(20) = {2+0.1, 4.5-0.1, 0, mesh_size_air};



// points to lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 5};

Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};

Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 17};

// generate surface
Coherence;
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8, 9, 10, 11, 12};
Curve Loop(3) = {13, 14, 15, 16};
Curve Loop(4) = {17, 18, 19, 20};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {1,2,3,4};
Plane Surface(4) = {4};

// physical surface
Physical Surface("group@s1@steel@demo_steel") = {1,2};
Color blue {Surface{1,2};}
Physical Surface("group@a1@other@air") = {3};
Color Red {Surface{3};}
Physical Surface("group@c1@copper@cu") = {4};
Color Orange {Surface{4};}

// boundary condition examples
Physical Line("boundary@dirichlet@0") = {1,2,4,10,11,12};
// Physical Line("boundary@dirichlet@1e-2") = {2};
// Physical Line("boundary@symmetry@odd1") = {2};
// Physical Line("boundary@symmetry@odd2") = {4};
// Physical Line("boundary@dirichlet@0") = {1,2,3,4};

Mesh 2;
Save "test-T_A_Lipo.msh";

//+gmsh spheres.geo -3 -format msh2
SetFactory("OpenCASCADE");

Cylinder(1) = {0, 0, 0, 0, 0, 1, 3, 2*Pi};
Cylinder(2) = {-1, 0, 0, 0, 0, 1, 0.5, 2*Pi};
Cylinder(3) = {1, 0, 0, 0, 0, 1, 0.5, 2*Pi};
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Volume{3}; }

//+
// Physical Surface("Top", 17) = {11, 5, 8};
// //+
// Physical Surface("Side", 18) = {10};
// //+
// Physical Surface("Bottom", 19) = {12, 6, 9};

Physical Surface("Boundaries", 17) = {12, 6, 9, 10, 11, 5, 8};
Physical Volume("Domain", 21) = {1};
Physical Volume("Left", 22) = {2};
Physical Volume("Right", 23) = {3};

//+
Transfinite Curve {4, 6, 9, 7, 4} = 15 Using Progression 1;

//+ gmsh spheres.geo -3 -format msh2
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 5, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(2) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(3) = {0, 0, 0, 2, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(4) = {0, 0, 0, 3, -Pi/2, Pi/2, 2*Pi};

// BooleanDifference(5)={Volume{4};}{Volume{3};};
// BooleanDifference(6)={Volume{3}; Delete;}{Volume{2};};
// BooleanDifference(7)={Volume{1}; Delete;}{Volume{4}; Delete;};

Physical Volume("Domain", 2) = {1, 2, 3, 4};

// Physical Surface("Buondary", 5) = {5};
// Physical Volume("Void", 6) = {1, 4, 2};
// Physical Volume("Outer", 7) = {4};
// Physical Volume("Inner", 8) = {2};//+
Physical Surface("Side", 4) = {1};

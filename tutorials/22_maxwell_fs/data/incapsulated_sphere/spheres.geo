
SetFactory("OpenCASCADE");
Sphere(1) = {-0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
Sphere(2) = {-0, 0, 0, 0.1, -Pi/2, Pi/2, 2*Pi};
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; }

Transfinite Curve {5} = 7 Using Progression 1;

Physical Volume("Big", 3) = {1};
Physical Volume("Small", 4) = {2};
Physical Surface("Boundary", 10) = {3};



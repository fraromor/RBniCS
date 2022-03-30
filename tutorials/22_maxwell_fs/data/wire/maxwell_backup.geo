SetFactory("OpenCASCADE");

// domain
Cylinder(1) = {0, 0, 0, 0, 0, 1, 5, 2*Pi};

// iron cylinder
Cylinder(2) = {0, 0.0, 0, 0, 0, 1, 1.2, 2*Pi};
Cylinder(3) = {0, 0.0, 0, 0, 0, 1, 1, 2*Pi};

// copper cylinders
n = 10;
c_1 = 0.8;
c_2 = 1.4;
r = 0.1;

For i In {1:n}
  angles_N~{i} = i*2*Pi/n;
  angles_S~{i} = (i + 0.5)*2*Pi/n;
EndFor

Cylinder(4) = {c_1*Cos(angles_N_1), c_1*Sin(angles_N_1), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(5) = {c_1*Cos(angles_N_2), c_1*Sin(angles_N_2), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(6) = {c_1*Cos(angles_N_3), c_1*Sin(angles_N_3), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(7) = {c_1*Cos(angles_N_4), c_1*Sin(angles_N_4), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(8) = {c_1*Cos(angles_N_5), c_1*Sin(angles_N_5), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(9) = {c_1*Cos(angles_N_6), c_1*Sin(angles_N_6), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(10) = {c_1*Cos(angles_N_7), c_1*Sin(angles_N_7), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(11) = {c_1*Cos(angles_N_8), c_1*Sin(angles_N_8), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(12) = {c_1*Cos(angles_N_9), c_1*Sin(angles_N_9), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(13) = {c_1*Cos(angles_N_10), c_1*Sin(angles_N_10), 0.0, 0, 0, 1, r, 2*Pi};

Cylinder(14) = {c_2*Cos(angles_S_1), c_2*Sin(angles_S_1), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(15) = {c_2*Cos(angles_S_2), c_2*Sin(angles_S_2), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(16) = {c_2*Cos(angles_S_3), c_2*Sin(angles_S_3), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(17) = {c_2*Cos(angles_S_4), c_2*Sin(angles_S_4), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(18) = {c_2*Cos(angles_S_5), c_2*Sin(angles_S_5), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(19) = {c_2*Cos(angles_S_6), c_2*Sin(angles_S_6), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(20) = {c_2*Cos(angles_S_7), c_2*Sin(angles_S_7), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(21) = {c_2*Cos(angles_S_8), c_2*Sin(angles_S_8), 0.0, 0, 0, 1, r,
2*Pi};
Cylinder(22) = {c_2*Cos(angles_S_9), c_2*Sin(angles_S_9), 0.0, 0, 0, 1, r, 2*Pi};
Cylinder(23) = {c_2*Cos(angles_S_10), c_2*Sin(angles_S_10), 0.0, 0, 0, 1, r,
  2*Pi};

BooleanDifference{Volume{2}; Delete;} {Volume{3}; Delete;}

BooleanDifference{Volume{1}; Delete;} { Volume{7}; Volume{6}; Volume{5};
Volume{4}; Volume{13}; Volume{12}; Volume{11}; Volume{10}; Volume{9}; Volume{8};
Volume{17}; Volume{16}; Volume{15}; Volume{14}; Volume{23}; Volume{22};
Volume{21}; Volume{20}; Volume{19}; Volume{18}; Volume{2};}

Physical Volume(70) = {24, 25};
Physical Volume(71) = {2};

Physical Volume(72) = {4:13};
Physical Volume(73) = {14:23};


SetFactory("OpenCASCADE");

// domain
Disk(1) = {0, 0, 0, 5, 5};

// iron Disk
Disk(2) = {0, 0.0, 0, 1.2, 1.2};
Disk(3) = {0, 0.0, 0, 1, 1};

// copper Disks
n = 10;
c_1 = 0.8;
c_2 = 1.4;
r = 0.1;

For i In {1:n}
  angles_N~{i} = i*2*Pi/n;
  angles_S~{i} = (i + 0.5)*2*Pi/n;
EndFor

Disk(4) = {c_1*Cos(angles_N_1), c_1*Sin(angles_N_1), 0.0, r, r};
Disk(5) = {c_1*Cos(angles_N_2), c_1*Sin(angles_N_2), 0.0, r, r};
Disk(6) = {c_1*Cos(angles_N_3), c_1*Sin(angles_N_3), 0.0, r, r};
Disk(7) = {c_1*Cos(angles_N_4), c_1*Sin(angles_N_4), 0.0, r, r};
Disk(8) = {c_1*Cos(angles_N_5), c_1*Sin(angles_N_5), 0.0, r, r};
Disk(9) = {c_1*Cos(angles_N_6), c_1*Sin(angles_N_6), 0.0, r, r};
Disk(10) = {c_1*Cos(angles_N_7), c_1*Sin(angles_N_7), 0.0, r, r};
Disk(11) = {c_1*Cos(angles_N_8), c_1*Sin(angles_N_8), 0.0, r, r};
Disk(12) = {c_1*Cos(angles_N_9), c_1*Sin(angles_N_9), 0.0, r, r};
Disk(13) = {c_1*Cos(angles_N_10), c_1*Sin(angles_N_10), 0.0, r, r};

Disk(14) = {c_2*Cos(angles_S_1), c_2*Sin(angles_S_1), 0.0, r, r};
Disk(15) = {c_2*Cos(angles_S_2), c_2*Sin(angles_S_2), 0.0, r, r};
Disk(16) = {c_2*Cos(angles_S_3), c_2*Sin(angles_S_3), 0.0, r, r};
Disk(17) = {c_2*Cos(angles_S_4), c_2*Sin(angles_S_4), 0.0, r, r};
Disk(18) = {c_2*Cos(angles_S_5), c_2*Sin(angles_S_5), 0.0, r, r};
Disk(19) = {c_2*Cos(angles_S_6), c_2*Sin(angles_S_6), 0.0, r, r};
Disk(20) = {c_2*Cos(angles_S_7), c_2*Sin(angles_S_7), 0.0, r, r};
Disk(21) = {c_2*Cos(angles_S_8), c_2*Sin(angles_S_8), 0.0, r, r};
Disk(22) = {c_2*Cos(angles_S_9), c_2*Sin(angles_S_9), 0.0, r, r};
Disk(23) = {c_2*Cos(angles_S_10), c_2*Sin(angles_S_10), 0.0, r, r};

BooleanDifference{Surface{2}; Delete;}{Surface{3}; Delete;}// Surface 2

BooleanDifference{ Surface{1}; Delete; }{ Surface{16}; Surface{15}; Surface{14};
Surface{23}; Surface{22}; Surface{21}; Surface{20}; Surface{19}; Surface{18};
Surface{17}; Surface{7}; Surface{6}; Surface{5}; Surface{4}; Surface{13};
Surface{12}; Surface{11}; Surface{10}; Surface{9}; Surface{8}; Surface{2};}// 24 25
//+

//+
Extrude {0, 0, 5} {
  Surface{24}; Surface{16}; Surface{17}; Surface{2}; Surface{18}; Surface{8}; Surface{7}; Surface{6}; Surface{15}; Surface{14}; Surface{5}; Surface{4}; Surface{23}; Surface{13}; Surface{22}; Surface{12}; Surface{11}; Surface{10}; Surface{9}; Surface{19}; Surface{20}; Surface{21}; Surface{25};
}

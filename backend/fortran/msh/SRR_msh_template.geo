// Template mesh geometry file for a split ring resonator,
// with its open side at the top.

d = 1; // grating period

d_in_nm  = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
lvert_nm = 0;
lhori_nm = 0;
width_nm = 0;

lvert = lvert_nm/d_in_nm*d;
lhori = lhori_nm/d_in_nm*d;
w     = width_nm/d_in_nm*d;

lc = 0;
lc2 = lc/1; // on SRR edgesq
lc3 = lc/1; //in centre of SRR

hy = dy; // Thickness: Squre profile => hy=d
hx = 0.;

Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0, lc};

// Vertices of the triangles

Point(5) = {-lhori/2+d/2,-lvert/2-d/2,0,lc2};
Point(6) = {-lhori/2+d/2,lvert/2-d/2,0,lc2};
Point(7) = {-lhori/2+d/2+w,lvert/2-d/2,0,lc2};
Point(8) = {-lhori/2+w+d/2,-lvert/2+w-d/2,0,lc2};
Point(9) = {lhori/2-w+d/2,-lvert/2+w-d/2,0,lc2};
Point(10) = {lhori/2-w+d/2,lvert/2-d/2,0,lc2};
Point(11) = {lhori/2+d/2,lvert/2-d/2,0,lc2};
Point(12) = {lhori/2+d/2,-lvert/2-d/2,0,lc2};
Point(13) = {-lhori/2+w+d/2,-lvert/2-d/2,0,lc2};
Point(14) = {lhori/2-w+d/2,-lvert/2-d/2,0,lc2};
Point(15) = {-lhori/2+d/2+w,-hy/2,0,lc2};
Point(16) = {lhori/2+d/2-w,-hy/2,0,lc2};
Point(17) = {d/2,lvert/2-d/2,0,lc2};
Point(18) = {d/2,-lvert/2-d/2+w,0,lc2};
Point(19) = {d/2,-hy/2,0,lc};


Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {6, 7};
Line(6) = {7, 15};
Line(7) = {15, 8};
Line(9) = {9, 16};
Line(10) = {16, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 14};
Line(14) = {14, 13};
Line(15) = {13, 5};
Line(16) = {5, 6};
Line(17) = {7, 17};
Line(18) = {17, 10};
Line(19) = {15, 19};
Line(20) = {19, 16};
Line(21) = {18, 19};
Line(22) = {19, 17};
Line(23) = {8, 13};
Line(24) = {9, 14};
Line(25) = {12, 3};
Line(26) = {6, 1};
Line(27) = {2, 5};
Line(28) = {8, 18};
Line(29) = {18, 9};
Line(30) = {4, 11};
Line Loop(31) = {26, 1, 30, -11, -18, -17, -5};
Plane Surface(32) = {31};
Line Loop(33) = {30, 12, 25, -2};
Plane Surface(34) = {33};
Line Loop(35) = {25, 3, 27, -15, -14, -13};
Plane Surface(36) = {35};
Line Loop(37) = {27, 16, 26, -4};
Plane Surface(38) = {37};
Line Loop(39) = {16, 5, 6, 7, 23, 15};
Plane Surface(40) = {39};
Line Loop(41) = {23, -14, -24, -29, -28};
Plane Surface(42) = {41};
Line Loop(43) = {24, -13, -12, -11, -10, -9};
Plane Surface(44) = {43};
Line Loop(45) = {29, 9, -20, -21};
Plane Surface(46) = {45};
Line Loop(47) = {21, -19, 7, 28};
Plane Surface(48) = {47};
Line Loop(49) = {6, 19, 22, -17};
Plane Surface(50) = {49};
Line Loop(51) = {22, 18, -10, -20};
Plane Surface(52) = {51};

Physical Line(53) = {4};
Physical Line(54) = {1};
Physical Line(56) = {3};
Physical Line(57) = {2};

Physical Surface(1) = {38, 32, 50, 52, 48, 46, 34, 36};
Physical Surface(2) = {40, 42, 44};

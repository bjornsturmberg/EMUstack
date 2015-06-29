// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0;
radius1 = (a1/(2*d_in_nm))*d;
ellipticity = 0;
square = 0;
strip = 0;
st = strip/d_in_nm;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Vertices of the triangles

Point(5) = {-hx+d/2, -hy/2, 0,lc3};

Point(6) = {-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(7) = {-hx+d/2-radius1, -hy/2, 0, lc2};
Point(8) = {-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(9) = {-hx+d/2+radius1, -hy/2, 0, lc2};

Point(10) = {-hx+d/2, 0, 0, lc};
Point(11) = {0,-hy/2, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc};
Point(13) = {d, -hy/2, 0, lc};
Point(14) = {0,-hy/2+st/2, 0, lc3};
Point(15) = {0,-hy/2-st/2, 0, lc3};
Point(16) = {d, -hy/2+st/2, 0, lc3};
Point(17) = {d, -hy/2-st/2, 0, lc3};
Point(18) = {-hx+d/2, -hy/2+st/2, 0, lc3};
Point(19) = {-hx+d/2, -hy/2-st/2, 0, lc3};

Line(1) = {1,10};
Line(2) = {10,4};
Line(3) = {2,12};
Line(4) = {12,3};
Line(21) = {1, 14};
Line(22) = {14, 11};
Line(23) = {11, 15};
Line(24) = {15, 2};
Line(25) = {11, 7};
Line(26) = {7, 5};
Line(27) = {5, 9};
Line(28) = {9, 13};
Line(29) = {4, 16};
Line(30) = {16, 13};
Line(31) = {13, 17};
Line(32) = {17, 3};
Line(33) = {14, 18};
Line(34) = {18, 16};
Line(35) = {15, 19};
Line(36) = {19, 17};
Line(37) = {10, 18};
Line(38) = {18, 6};
Line(39) = {6, 5};
Line(40) = {5, 8};
Line(41) = {8, 19};
Line(42) = {19, 12};


If(square == 0)
Ellipsis(17) = {9,5,6,6};
Ellipsis(18) = {6,5,7,7};
Ellipsis(19) = {7,5,8,8};
Ellipsis(20) = {8,5,9,9};

Line Loop(43) = {21, 33, -37, -1};
Plane Surface(44) = {43};
Line Loop(45) = {37, 34, -29, -2};
Plane Surface(46) = {45};
Line Loop(47) = {24, 3, -42, -35};
Plane Surface(48) = {47};
Line Loop(49) = {42, 4, -32, -36};
Plane Surface(50) = {49};
Line Loop(51) = {22, 25, -18, -38, -33};
Plane Surface(52) = {51};
Line Loop(53) = {23, 35, -41, -19, -25};
Plane Surface(54) = {53};
Line Loop(55) = {38, -17, 28, -30, -34};
Plane Surface(56) = {55};
Line Loop(57) = {28, 31, -36, -41, 20};
Plane Surface(58) = {57};
Line Loop(59) = {19, -40, -26};
Plane Surface(60) = {59};
Line Loop(61) = {18, 26, -39};
Plane Surface(62) = {61};
Line Loop(63) = {39, 27, 17};
Plane Surface(64) = {63};
Line Loop(65) = {20, -27, 40};
Plane Surface(66) = {65};

Physical Line(67) = {21, 22, 23, 24};
Physical Line(68) = {1, 2};
Physical Line(69) = {29, 30, 31, 32};
Physical Line(70) = {4, 3};

Physical Surface(1) = {44, 46, 48, 50};
Physical Surface(2) = {62, 64, 66, 60};
Physical Surface(3) = {52, 56, 58, 54};
EndIf

If(square == 1)
// square
Point(150) = {-hx+d/2+radius1, -hy/2+radius1, 0,lc3};
Point(151) = {-hx+d/2-radius1, -hy/2+radius1, 0,lc3};
Point(152) = {-hx+d/2+radius1, -hy/2-radius1, 0,lc3};
Point(153) = {-hx+d/2-radius1, -hy/2-radius1, 0,lc3};

Line(51) = {7, 151};
Line(52) = {151, 6};
Line(53) = {6, 150};
Line(54) = {150, 9};
Line(55) = {7, 153};
Line(56) = {153, 8};
Line(57) = {152, 8};
Line(58) = {152, 9};

Line Loop(43) = {21, 33, -37, -1};
Plane Surface(44) = {43};
Line Loop(45) = {37, 34, -29, -2};
Plane Surface(46) = {45};
Line Loop(47) = {24, 3, -42, -35};
Plane Surface(48) = {47};
Line Loop(49) = {42, 4, -32, -36};
Plane Surface(50) = {49};
Line Loop(59) = {22, 25, 51, 52, -38, -33};
Plane Surface(60) = {59};
Line Loop(61) = {38, 53, 54, 28, -30, -34};
Plane Surface(62) = {61};
Line Loop(63) = {58, 28, 31, -36, -41, -57};
Plane Surface(64) = {63};
Line Loop(65) = {41, -35, -23, 25, 55, 56};
Plane Surface(66) = {65};
Line Loop(67) = {51, 52, 39, -26};
Plane Surface(68) = {67};
Line Loop(69) = {39, 27, -54, -53};
Plane Surface(70) = {69};
Line Loop(71) = {40, -57, 58, -27};
Plane Surface(72) = {71};
Line Loop(73) = {26, 40, -56, -55};
Plane Surface(74) = {73};

Physical Surface(1) = {44, 46, 48, 50};
Physical Surface(2) = {68, 70, 72, 74};
Physical Surface(3) = {60, 62, 64, 66};

EndIf




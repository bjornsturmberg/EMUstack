// Template mesh geometry file for a hexagonal unitcell.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0;
radius1 = (a1/(2*d_in_nm))*d;
radius2 = (a1/(2*d_in_nm))*d;
radius3 = (a1/(2*d_in_nm))*d;
radius4 = (a1/(2*d_in_nm))*d;
radius5 = (a1/(2*d_in_nm))*d;
ellipticity = 0;
square = 0;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres
lc4 = lc/1; // centres of top and bottom

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0, lc};

// Vertices of the triangles

Point(5) = {-hx+d/2, -hy/2, 0,lc3};

Point(6) = {-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(7) = {-hx+d/2-radius1, -hy/2, 0, lc2};
Point(8) = {-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(9) = {-hx+d/2+radius1, -hy/2, 0, lc2};

Point(10) = {-hx+d/2, 0, 0, lc4};
Point(11) = {0,-hy/2, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc4};
Point(13) = {d, -hy/2, 0, lc};

Point(14) = {0, -radius2, 0, lc2};
Point(15) = {radius2, 0, 0, lc2};
Point(16) = {-hx+radius3, -hy, 0, lc2};
Point(17) = {-hx, -hy+radius3, 0, lc2};
Point(18) = {-hx+d-radius4, -hy, 0, lc2};
Point(19) = {-hx+d, -hy+radius4, 0, lc2};
Point(20) = {d-radius5, 0, 0, lc2};
Point(21) = {d, -radius5, 0, lc2};


Line(1) = {15,10};
Line(2) = {20,4};
Line(3) = {2,16};
Line(4) = {18,3};
Line(5) = {14,11};
Line(6) = {17,2};
Line(7) = {4,21};
Line(8) = {19,3};
Line(9) = {11,7};
Line(10) = {7,5};
Line(11) = {5,9};
Line(12) = {9,13};
Line(13) = {10,6};
Line(14) = {6,5};
Line(15) = {5,8};
Line(16) = {8,12};
Line(49) = {1, 14};
Line(50) = {1, 15};
Line(51) = {10, 20};
Line(52) = {21, 13};
Line(53) = {13, 19};
Line(54) = {18, 12};
Line(55) = {12, 16};
Line(56) = {17, 11};

If(square == 0)
Ellipsis(17) = {9,5,6,6};
Ellipsis(18) = {6,5,7,7};
Ellipsis(19) = {7,5,8,8};
Ellipsis(20) = {8,5,9,9};
Ellipse(57) = {14, 1, 15, 15};
Ellipse(58) = {20, 4, 21, 21};
Ellipse(59) = {19, 3, 18, 18};
Ellipse(60) = {16, 2, 17, 17};

//Physical Line(61) = {50, 1, 51, 2};
//Physical Line(62) = {7, 52, 53, 8};
//Physical Line(63) = {4, 54, 55, 3};
//Physical Line(64) = {6, 56, 5, 49};

Line Loop(65) = {1, 13, 18, -9, -5, 57};
Plane Surface(66) = {65};
Line Loop(67) = {50, -57, -49};
Plane Surface(68) = {67};
Line Loop(69) = {51, 58, 52, -12, 17, -13};
Plane Surface(70) = {69};
Line Loop(71) = {7, -58, 2};
Plane Surface(72) = {71};
Line Loop(73) = {53, 59, 54, -16, 20, 12};
Plane Surface(74) = {73};
Line Loop(75) = {59, 4, -8};
Plane Surface(76) = {75};
Line Loop(77) = {20, -11, 15};
Plane Surface(78) = {77};
Line Loop(79) = {17, 14, 11};
Plane Surface(80) = {79};
Line Loop(81) = {18, 10, -14};
Plane Surface(82) = {81};
Line Loop(83) = {19, -15, -10};
Plane Surface(84) = {83};
Line Loop(85) = {16, 55, 60, 56, 9, 19};
Plane Surface(86) = {85};
Line Loop(87) = {60, 6, 3};
Plane Surface(88) = {87};

Physical Surface(1) = {66, 70, 74, 86};
Physical Surface(2) = {82, 80, 78, 84, 68, 72, 76, 88};

EndIf


If(square == 1)
// square
Point(150) = {-hx+d/2+radius1, -hy/2+radius1, 0,lc3};
Point(151) = {-hx+d/2-radius1, -hy/2+radius1, 0,lc3};
Point(152) = {-hx+d/2+radius1, -hy/2-radius1, 0,lc3};
Point(153) = {-hx+d/2-radius1, -hy/2-radius1, 0,lc3};
Point(154) = {radius2, -radius2, 0, lc3};
Point(156) = {-hx+radius3, -hy+radius3, 0, lc3};
Point(158) = {-hx+d-radius4, -hy+radius4, 0, lc3};
Point(159) = {d-radius5, -radius5, 0, lc3};


Line(17) = {151, 6};
Line(18) = {6, 150};
Line(19) = {150, 9};
Line(20) = {9, 152};
Line(21) = {152, 8};
Line(22) = {8, 153};
Line(23) = {153, 7};
Line(24) = {7, 151};
Line(41) = {14, 154};
Line(42) = {154, 15};
Line(43) = {20, 159};
Line(44) = {159, 21};
Line(45) = {19, 158};
Line(46) = {158, 18};
Line(47) = {17, 156};
Line(48) = {156, 16};

//Physical Line(61) = {50, 1, 51, 2};
//Physical Line(62) = {7, 52, 53, 8};
//Physical Line(63) = {4, 54, 55, 3};
//Physical Line(64) = {6, 56, 5, 49};

Line Loop(65) = {49, 41, 42, -50};
Plane Surface(66) = {65};
Line Loop(67) = {5, 9, 24, 17, -13, -1, -42, -41};
Plane Surface(68) = {67};
Line Loop(69) = {56, 9, -23, -22, 16, 55, -48, -47};
Plane Surface(70) = {69};
Line Loop(71) = {48, -3, -6, 47};
Plane Surface(72) = {71};
Line Loop(73) = {24, 17, 14, -10};
Plane Surface(74) = {73};
Line Loop(75) = {10, 15, 22, 23};
Plane Surface(76) = {75};
Line Loop(77) = {15, -21, -20, -11};
Plane Surface(78) = {77};
Line Loop(79) = {19, -11, -14, 18};
Plane Surface(80) = {79};
Line Loop(81) = {13, 18, 19, 12, -52, -44, -43, -51};
Plane Surface(82) = {81};
Line Loop(83) = {44, -7, -2, 43};
Plane Surface(84) = {83};
Line Loop(85) = {12, 53, 45, 46, 54, -16, -21, -20};
Plane Surface(86) = {85};
Line Loop(87) = {45, 46, 4, -8};
Plane Surface(88) = {87};
Physical Surface(1) = {68, 82, 86, 70};
Physical Surface(2) = {74, 80, 78, 76, 66, 84, 88, 72};

EndIf

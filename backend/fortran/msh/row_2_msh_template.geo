// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 10;
radius1 = (a1/(2*d_in_nm))*d;
sep = 20;
sep_norm = sep/d_in_nm;
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
Point(4) = {d, 0, 0,lc};

// Vertices of the triangles

Point(5) = {-sep_norm/2-hx+d/2, -hy/2, 0,lc3};
Point(6) = {-sep_norm/2-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(7) = {-sep_norm/2-hx+d/2-radius1, -hy/2, 0, lc2};
Point(8) = {-sep_norm/2-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(9) = {-sep_norm/2-hx+d/2+radius1, -hy/2, 0, lc2};

Point(10) = {-sep_norm/2-hx+d/2, 0, 0, lc4};
Point(11) = {0,-hy/2, 0, lc};
Point(12) = {-sep_norm/2-hx+d/2, -hy, 0, lc4};
Point(13) = {d, -hy/2, 0, lc};

Point(14) = {sep_norm/2-hx+d/2, -hy/2, 0,lc3};
Point(15) = {sep_norm/2-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(16) = {sep_norm/2-hx+d/2-radius1, -hy/2, 0, lc2};
Point(17) = {sep_norm/2-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(18) = {sep_norm/2-hx+d/2+radius1, -hy/2, 0, lc2};

Point(19) = {sep_norm/2-hx+d/2, 0, 0, lc4};
Point(20) = {sep_norm/2-hx+d/2, -hy, 0, lc4};

Line(1) = {1,10};
Line(2) = {10,19};
Line(3) = {2,12};
Line(4) = {12,20};
Line(5) = {1,11};
Line(6) = {11,2};
Line(7) = {4,13};
Line(8) = {13,3};
Line(10) = {7,5};
Line(11) = {5,9};
Line(14) = {6,5};
Line(15) = {5,8};
Line(57) = {19, 4};
Line(58) = {3, 20};

Line(21) = {11, 7};
Line(22) = {10, 6};
Line(23) = {8, 12};
Line(24) = {9, 16};
Line(25) = {18, 13};
Line(26) = {19, 15};
Line(27) = {17, 20};
Line(32) = {15, 14};
Line(33) = {14, 17};
Line(34) = {14, 18};
Line(35) = {14, 16};

Physical Line(59) = {1, 2, 57};
Physical Line(60) = {7, 8};
Physical Line(61) = {58, 4, 3};
Physical Line(62) = {6, 5};

//If(ellipticity == 0)
//Circle(17) = {9,5,6};
//Circle(18) = {6,5,7};
//Circle(19) = {7,5,8};
//Circle(20) = {8,5,9};
//EndIf
If(square == 0)
Ellipsis(17) = {9,5,6,6};
Ellipsis(18) = {6,5,7,7};
Ellipsis(19) = {7,5,8,8};
Ellipsis(20) = {8,5,9,9};
Ellipse(28) = {18, 14, 15, 15};
Ellipse(29) = {18, 14, 17, 17};
Ellipse(30) = {17, 14, 16, 16};
Ellipse(31) = {16, 14, 15, 15};


Line Loop(63) = {21, -18, -22, -1, 5};
Plane Surface(64) = {63};
Line Loop(65) = {22, -17, 24, 31, -26, -2};
Plane Surface(66) = {65};
Line Loop(67) = {57, 7, -25, 28, -26};
Plane Surface(68) = {67};
Line Loop(69) = {25, 8, 58, -27, -29};
Plane Surface(70) = {69};
Line Loop(71) = {27, -4, -23, 20, 24, -30};
Plane Surface(72) = {71};
Line Loop(73) = {23, -3, -6, 21, 19};
Plane Surface(74) = {73};
Line Loop(75) = {19, -15, -10};
Plane Surface(76) = {75};
Line Loop(77) = {18, 10, -14};
Plane Surface(78) = {77};
Line Loop(79) = {14, 11, 17};
Plane Surface(80) = {79};
Line Loop(81) = {20, -11, 15};
Plane Surface(82) = {81};
Line Loop(83) = {30, -35, 33};
Plane Surface(84) = {83};
Line Loop(85) = {35, 31, 32};
Plane Surface(86) = {85};
Line Loop(87) = {32, 34, 28};
Plane Surface(88) = {87};
Line Loop(89) = {34, 29, -33};
Plane Surface(90) = {89};
Physical Surface(1) = {64, 66, 68, 70, 72, 74};
Physical Surface(2) = {78, 80, 82, 76};
Physical Surface(3) = {86, 88, 90, 84};

EndIf
If(square == 1)
// square
Point(150) = {-sep_norm/2-hx+d/2+radius1, -hy/2+radius1, 0,lc3};
Point(151) = {-sep_norm/2-hx+d/2-radius1, -hy/2+radius1, 0,lc3};
Point(152) = {-sep_norm/2-hx+d/2+radius1, -hy/2-radius1, 0,lc3};
Point(153) = {-sep_norm/2-hx+d/2-radius1, -hy/2-radius1, 0,lc3};
Point(154) = {sep_norm/2-hx+d/2+radius1, -hy/2+radius1, 0,lc3};
Point(155) = {sep_norm/2-hx+d/2-radius1, -hy/2+radius1, 0,lc3};
Point(156) = {sep_norm/2-hx+d/2+radius1, -hy/2-radius1, 0,lc3};
Point(157) = {sep_norm/2-hx+d/2-radius1, -hy/2-radius1, 0,lc3};


Line(21) = {152, 8};
Line(22) = {8, 153};
Line(23) = {153, 7};
Line(24) = {7, 151};
Line(41) = {151, 6};
Line(42) = {6, 150};
Line(43) = {150, 9};
Line(44) = {151, 7};
Line(45) = {7, 153};
Line(46) = {153, 8};
Line(47) = {8, 152};
Line(48) = {152, 9};
Line(49) = {155, 15};
Line(50) = {15, 154};
Line(51) = {154, 18};
Line(52) = {16, 155};
Line(53) = {157, 17};
Line(54) = {17, 156};
Line(55) = {156, 18};
Line(56) = {16, 157};


Line Loop(63) = {1, 22, -41, 44, -21, -5};
Plane Surface(64) = {63};
Line Loop(65) = {22, 42, 43, 24, 52, 49, -26, -2};
Plane Surface(66) = {65};
Line Loop(67) = {26, 50, 51, 25, -7, -57};
Plane Surface(68) = {67};
Line Loop(69) = {54, 55, 25, 8, 58, -27};
Plane Surface(70) = {69};
Line Loop(71) = {27, -4, -23, 47, 48, 24, 56, 53};
Plane Surface(72) = {71};
Line Loop(73) = {23, -3, -6, 21, 45, 46};
Plane Surface(74) = {73};
Line Loop(75) = {41, 14, -10, -44};
Plane Surface(76) = {75};
Line Loop(77) = {14, 11, -43, -42};
Plane Surface(78) = {77};
Line Loop(79) = {15, 47, 48, -11};
Plane Surface(80) = {79};
Line Loop(81) = {15, -46, -45, 10};
Plane Surface(82) = {81};
Line Loop(83) = {56, 53, -33, 35};
Line Loop(84) = {52, 49, 32, 35};
Plane Surface(85) = {84};
Line Loop(86) = {32, 34, -51, -50};
Plane Surface(87) = {86};
Line Loop(88) = {55, -34, 33, 54};
Plane Surface(89) = {88};
Plane Surface(90) = {83};
Physical Surface(1) = {64, 66, 68, 70, 72, 74};
Physical Surface(2) = {76, 78, 82, 80};
Physical Surface(3) = {85, 87, 89, 90};
EndIf

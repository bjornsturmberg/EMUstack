// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0; // outside radius
radius1 = (a1/(2*d_in_nm))*d;
a2 = 0; // inside radius
radius2 = (a2/(2*d_in_nm))*d;
ellipticity = 0;
square = 0;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres

xshift_nm = 0;
xshift = xshift_nm/d_in_nm;

hy = dy; // Thickness: Squre profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Vertices of the triangles

Point(5) = {xshift-hx+d/2, -hy/2, 0,lc3};

Point(6) = {xshift-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(7) = {xshift-hx+d/2-radius1, -hy/2, 0, lc2};
Point(8) = {xshift-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(9) = {xshift-hx+d/2+radius1, -hy/2, 0, lc2};

Point(10) = {xshift-hx+d/2, 0, 0, lc};
Point(11) = {0,-hy/2, 0, lc};
Point(12) = {xshift-hx+d/2, -hy, 0, lc};
Point(13) = {d, -hy/2, 0, lc};

Point(14) = {xshift-hx+d/2, -hy/2+(radius2-ellipticity*radius2), 0, lc3};
Point(15) = {xshift-hx+d/2-radius2, -hy/2, 0, lc3};
Point(16) = {xshift-hx+d/2, -hy/2-(radius2-ellipticity*radius2), 0, lc3};
Point(17) = {xshift-hx+d/2+radius2, -hy/2, 0, lc3};

Line(1) = {1, 10};
Line(2) = {10, 4};
Line(3) = {4, 13};
Line(4) = {13, 3};
Line(5) = {3, 12};
Line(6) = {12, 2};
Line(7) = {2, 11};
Line(8) = {11, 1};
Line(9) = {11, 7};
Line(10) = {7, 15};
Line(11) = {15, 5};
Line(12) = {10, 6};
Line(13) = {6, 14};
Line(14) = {14, 5};
Line(15) = {12, 8};
Line(16) = {8, 16};
Line(17) = {16, 5};
Line(18) = {13, 9};
Line(19) = {9, 17};
Line(20) = {17, 5};


If(square == 0)
Ellipse(21) = {7, 5, 6, 6};
Ellipse(22) = {6, 5, 9, 9};
Ellipse(23) = {9, 5, 8, 8};
Ellipse(24) = {8, 5, 7, 7};
Ellipse(25) = {15, 5, 14, 14};
Ellipse(26) = {14, 5, 17, 17};
Ellipse(27) = {17, 5, 16, 16};
Ellipse(28) = {16, 5, 15, 15};

Line Loop(29) = {7, 9, -24, -15, 6};
Plane Surface(30) = {29};
Line Loop(31) = {9, 21, -12, -1, -8};
Plane Surface(32) = {31};
Line Loop(33) = {12, 22, -18, -3, -2};
Plane Surface(34) = {33};
Line Loop(35) = {18, 23, -15, -5, -4};
Plane Surface(36) = {35};
Line Loop(37) = {25, 14, -11};
Plane Surface(38) = {37};
Line Loop(39) = {14, -20, -26};
Plane Surface(40) = {39};
Line Loop(41) = {20, -17, -27};
Plane Surface(42) = {41};
Line Loop(43) = {17, -11, -28};
Plane Surface(44) = {43};
Line Loop(45) = {21, 13, -25, -10};
Plane Surface(46) = {45};
Line Loop(47) = {13, 26, -19, -22};
Plane Surface(48) = {47};
Line Loop(49) = {27, -16, -23, 19};
Plane Surface(50) = {49};
Line Loop(51) = {16, 28, -10, -24};
Plane Surface(52) = {51};

Physical Line(53) = {1, 2};
Physical Line(54) = {3, 4};
Physical Line(55) = {5, 6};
Physical Line(56) = {7, 8};

Physical Surface(1) = {32, 34, 36, 30};
Physical Surface(2) = {46, 48, 50, 52};
Physical Surface(3) = {44, 38, 40, 42};



EndIf
If(square == 1)
// square
Point(150) = {-hx+d/2+radius1, -hy/2+radius1, 0,lc2};
Point(151) = {-hx+d/2-radius1, -hy/2+radius1, 0,lc2};
Point(152) = {-hx+d/2+radius1, -hy/2-radius1, 0,lc2};
Point(153) = {-hx+d/2-radius1, -hy/2-radius1, 0,lc2};
Point(154) = {-hx+d/2+radius2, -hy/2+radius2, 0,lc2};
Point(155) = {-hx+d/2-radius2, -hy/2+radius2, 0,lc2};
Point(156) = {-hx+d/2+radius2, -hy/2-radius2, 0,lc2};
Point(157) = {-hx+d/2-radius2, -hy/2-radius2, 0,lc2};

Line(21) = {151, 6};
Line(22) = {6, 150};
Line(23) = {150, 9};
Line(24) = {9, 152};
Line(25) = {152, 8};
Line(26) = {8, 153};
Line(27) = {153, 7};
Line(28) = {7, 151};
Line(29) = {155, 15};
Line(30) = {15, 157};
Line(31) = {157, 16};
Line(32) = {16, 156};
Line(33) = {156, 17};
Line(34) = {17, 154};
Line(35) = {154, 14};
Line(36) = {14, 155};

Line Loop(37) = {8, 1, 12, -21, -28, -9};
Plane Surface(38) = {37};
Line Loop(39) = {12, 22, 23, -18, -3, -2};
Plane Surface(40) = {39};
Line Loop(41) = {18, 24, 25, -15, -5, -4};
Plane Surface(42) = {41};
Line Loop(43) = {15, 26, 27, -9, -7, -6};
Plane Surface(44) = {43};
Line Loop(45) = {29, 11, -14, 36};
Plane Surface(46) = {45};
Line Loop(47) = {14, -20, 34, 35};
Plane Surface(48) = {47};
Line Loop(49) = {20, -17, 32, 33};
Plane Surface(50) = {49};
Line Loop(51) = {17, -11, 30, 31};
Plane Surface(52) = {51};
Line Loop(53) = {28, 21, 13, 36, 29, -10};
Plane Surface(54) = {53};
Line Loop(55) = {13, -35, -34, -19, -23, -22};
Plane Surface(56) = {55};
Line Loop(57) = {19, -33, -32, -16, -25, -24};
Plane Surface(58) = {57};
Line Loop(59) = {16, -31, -30, -10, -27, -26};
Plane Surface(60) = {59};

Physical Line(53) = {1, 2};
Physical Line(54) = {3, 4};
Physical Line(55) = {5, 6};
Physical Line(56) = {7, 8};


Physical Surface(1) = {38, 40, 42, 44, 46, 52, 50, 48};
Physical Surface(2) = {54, 56, 58, 60};

EndIf

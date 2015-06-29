d = 1; // grating period
ff = 0;
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0;
a2 = 0;
radius1 = (a1/(2*d_in_nm))*d;
radius2 = (a2/(2*d_in_nm))*d;
//radius2 = ((ff*(d)^2)/3.14159265 - (radius1^2))^0.5;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres
lc4 = lc/1; // cylinder2 centres
posx = 0;//
posy = 0;//

hy = dy; // Thickness: Squre profile => hy=d
hx = 0.;

// 2*2 supercell outline

Point(1) = {0, 0, 0, lc4};
Point(2) = {-hx, -hy, 0, lc4};
Point(3) = {-hx+d, -hy, 0, lc4};
Point(4) = {d, 0, 0,lc4};
Point(5) = {-hx+d/2+posx, -hy/2+posy, 0,lc3};
Point(10) = {-hx+d/2+posx, 0, 0, lc};
Point(11) = {0,-hy/2+posy, 0, lc};
Point(12) = {-hx+d/2+posx, -hy, 0, lc};
Point(13) = {d, -hy/2+posy, 0, lc};

Point(14) = {-hx+radius2, 0, 0, lc2};
Point(15) = {d-radius2, 0, 0, lc2};
Point(16) = {0,-radius2, 0, lc2};
Point(18) = {-hx+d/2+posx, -hy/2+radius1+posy, 0, lc2};
Point(20) = {d,-radius2, 0, lc2};
Point(21) = {-hx+d/2-radius1+posx, -hy/2+posy, 0,lc2};
Point(22) = {-hx+d/2+radius1+posx, -hy/2+posy, 0,lc2};
Point(23) = {-hx, -hy+radius2, 0, lc2};
Point(25) = {-hx+d/2+posx, -hy/2-radius1+posy, 0, lc2};
Point(27) = {d,-hy+radius2, 0, lc2};
Point(28) = {-hx+radius2, -hy, 0, lc2};
Point(29) = {d-radius2, -hy, 0, lc2};

Line(1) = {1, 14};
Line(2) = {14, 10};
Line(3) = {10, 15};
Line(4) = {15, 4};
Line(5) = {4, 20};
Line(6) = {20, 13};
Line(7) = {13, 27};
Line(8) = {27, 3};
Line(9) = {3, 29};
Line(10) = {29, 12};
Line(11) = {12, 28};
Line(12) = {28, 2};
Line(13) = {2, 23};
Line(14) = {23, 11};
Line(15) = {11, 16};
Line(16) = {16, 1};
Line(18) = {18, 10};
Line(19) = {18, 5};
Line(20) = {5, 25};
Line(21) = {25, 12};
Line(22) = {11, 21};
Line(23) = {21, 5};
Line(24) = {5, 22};
Line(25) = {22, 13};

Ellipsis(26) = {14, 1, 16, 16};
Ellipsis(27) = {15, 4, 20, 20};
Ellipsis(28) = {27, 3, 29, 29};
Ellipsis(29) = {28, 2, 23, 23};
Ellipsis(30) = {21, 5, 18, 18};
Ellipsis(31) = {18, 5, 22, 22};
Ellipsis(32) = {22, 5, 25, 25};
Ellipsis(33) = {25, 5, 21, 21};

Line Loop(77) = {1, 26, 16};
Plane Surface(78) = {77};
Line Loop(79) = {4, 5, -27};
Plane Surface(80) = {79};
Line Loop(81) = {28, -9, -8};
Plane Surface(82) = {81};
Line Loop(83) = {29, -13, -12};
Plane Surface(84) = {83};
Line Loop(85) = {30, 19, -23};
Plane Surface(86) = {85};
Line Loop(87) = {19, 24, -31};
Plane Surface(88) = {87};
Line Loop(89) = {24, 32, -20};
Plane Surface(90) = {89};
Line Loop(91) = {20, 33, 23};
Plane Surface(92) = {91};
Line Loop(93) = {2, -18, -30, -22, 15, -26};
Plane Surface(94) = {93};
Line Loop(95) = {3, 27, 6, -25, -31, 18};
Plane Surface(96) = {95};
Line Loop(97) = {25, 7, 28, 10, -21, -32};
Plane Surface(98) = {97};
Line Loop(99) = {21, 11, 29, 14, 22, -33};
Plane Surface(100) = {99};

Physical Line(101) = {1, 2, 3, 4};
Physical Line(102) = {5, 6, 7, 8};
Physical Line(103) = {9, 10, 11, 12};
Physical Line(104) = {13, 14, 15, 16};

Physical Surface(1) = {94, 96, 98, 100};
Physical Surface(2) = {86, 88, 90, 92};
Physical Surface(3) = {78, 80, 82, 84};

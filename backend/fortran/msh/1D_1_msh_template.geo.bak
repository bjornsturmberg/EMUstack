d = 1; // unit cell period
d_in_nm = 0; // grating period
w1 = 0;
h_width1 = (w1/(2*d_in_nm))*d;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on grating surfaces

hy = d; // Thickness: Squre profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Edges of grating
Point(6) = {hy/2-h_width1, 0, 0, lc2};
Point(7) = {hy/2+h_width1, 0, 0, lc2};
Point(8) = {hy/2+h_width1, -hy, 0, lc2};
Point(9) = {hy/2-h_width1, -hy, 0, lc2};

//Centre of grating
Point(10) = {hy/2, 0, 0, lc};
Point(11) = {hy/2, -d, 0, lc};

Line(1) = {1, 6};
Line(2) = {6, 10};
Line(3) = {10, 7};
Line(4) = {7, 4};
Line(5) = {4, 3};
Line(6) = {3, 8};
Line(7) = {8, 11};
Line(8) = {11, 9};
Line(9) = {9, 2};
Line(10) = {2, 1};
Line(11) = {6, 9};
Line(12) = {10, 11};
Line(13) = {7, 8};

Line Loop(14) = {1, 11, 9, 10};
Plane Surface(15) = {14};
Line Loop(16) = {13, -6, -5, -4};
Plane Surface(17) = {16};
Line Loop(18) = {11, -8, -12, -2};
Plane Surface(19) = {18};
Line Loop(20) = {12, -7, -13, -3};
Plane Surface(21) = {20};


Physical Line(22) = {10};
Physical Line(23) = {1, 2, 3, 4};
Physical Line(24) = {5};
Physical Line(25) = {6, 7, 8, 9};


Physical Surface(1) = {15, 17};
Physical Surface(2) = {19, 21};

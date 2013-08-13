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
Point(6) = {0, -hy/2-h_width1, 0, lc2};
Point(7) = {0, -hy/2+h_width1, 0, lc2};
Point(8) = {d, -hy/2+h_width1, 0, lc2};
Point(9) = {d, -hy/2-h_width1, 0, lc2};

//Centre of grating
Point(10) = {0, -hy/2, 0, lc};
Point(11) = {d, -hy/2, 0, lc};

Line(1) = {1, 4};
Line(2) = {4, 8};
Line(3) = {8, 11};
Line(4) = {11, 9};
Line(5) = {9, 3};
Line(6) = {3, 2};
Line(7) = {2, 6};
Line(8) = {6, 10};
Line(9) = {10, 7};
Line(10) = {7, 1};
Line(11) = {7, 8};
Line(12) = {10, 11};
Line(13) = {6, 9};


Line Loop(14) = {10, 1, 2, -11};
Plane Surface(15) = {14};
Line Loop(16) = {11, 3, -12, 9};
Plane Surface(17) = {16};
Line Loop(18) = {8, 12, 4, -13};
Plane Surface(19) = {18};
Line Loop(20) = {7, 13, 5, 6};
Plane Surface(21) = {20};


Physical Line(22) = {1};
Physical Line(23) = {2, 3, 4, 5};
Physical Line(24) = {6};
Physical Line(25) = {7, 8, 9, 10};


Physical Surface(1) = {15, 21};
Physical Surface(2) = {17, 19};

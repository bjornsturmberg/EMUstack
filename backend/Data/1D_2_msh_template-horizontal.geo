d = 1; // grating period
d_in_nm = 0;
w1 = 0;
w2 = 0;
h_width1 = (w1/(2*d_in_nm))*d;
h_width2 = (w2/(2*d_in_nm))*d;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on grating_1 surfaces
lc3 = lc/1; // on grating_2 surfaces

hy = d; // Thickness: Squre profile => hy=d
hx = 0.;

posx1 = hy/4;// distance from bottom to bottom grating center
posx2 = hy/4;// distance from top to top grating center
posx3 = hy/2;// middle point betwen top & bottom grating centers


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

Point(5) = {-hx, -hy+posx3, 0,lc};
Point(6) = {d,   -hy+posx3, 0,lc};

// Edges of grating
Point(7) =  {0, -hy+posx1-h_width1, 0, lc2};
Point(8) =  {0, -hy+posx1+h_width1, 0, lc2};
Point(9) =  {d, -hy+posx1+h_width1, 0, lc2};
Point(10) = {d, -hy+posx1-h_width1, 0, lc2};

Point(11) = {0, 0-posx2-h_width2, 0, lc3};
Point(12) = {0, 0-posx2+h_width2, 0, lc3};
Point(13) = {d, 0-posx2+h_width2, 0, lc3};
Point(14) = {d, 0-posx2-h_width2, 0, lc3};



Line(1) = {1, 4};
Line(2) = {4, 13};
Line(3) = {13, 14};
Line(4) = {14, 6};
Line(5) = {6, 9};
Line(6) = {9, 10};
Line(7) = {10, 3};
Line(8) = {3, 2};
Line(9) = {2, 7};
Line(10) = {8, 7};
Line(11) = {8, 5};
Line(12) = {5, 11};
Line(13) = {11, 12};
Line(14) = {12, 1};
Line(15) = {12, 13};
Line(16) = {11, 14};
Line(17) = {5, 6};
Line(18) = {8, 9};
Line(19) = {7, 10};


Line Loop(20) = {1, 2, -15, 14};
Plane Surface(21) = {20};
Line Loop(22) = {15, 3, -16, 13};
Plane Surface(23) = {22};
Line Loop(24) = {16, 4, -17, 12};
Plane Surface(25) = {24};
Line Loop(26) = {11, 17, 5, -18};
Plane Surface(27) = {26};
Line Loop(28) = {10, 19, -6, -18};
Plane Surface(29) = {28};
Line Loop(30) = {19, 7, 8, 9};
Plane Surface(31) = {30};


Physical Line(32) = {1};
Physical Line(33) = {2, 3, 4, 5, 6, 7};
Physical Line(34) = {8};
Physical Line(35) = {9, 10, 11, 12, 13, 14};


Physical Surface(1) = {21, 25, 27, 31};
Physical Surface(2) = {23};
Physical Surface(3) = {29};

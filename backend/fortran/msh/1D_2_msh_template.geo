d = 1; // grating period
d_in_nm = 0;
w1 = 0;
w2 = 0;
h_width1 = (w1/(2*d_in_nm))*d;
h_width2 = (w2/(2*d_in_nm))*d;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on grating_1 surfaces
lc3 = lc/1; // on grating_1 center
lc4 = lc/1; // on grating_2 center

hy = d; // Thickness: Squre profile => hy=d
hx = 0.;

posx1 = hy/4;// distance from bottom to bottom grating center
posx2 = hy/4;// distance from top to top grating center
posx3 = hy/2;// middle point betwen top & bottom grating centers


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

Point(5) = {hx+posx3, hx, 0,lc};
Point(6) = {hx+posx3, -d,   0,lc};

// Edges of grating
Point(7) =  {posx3+posx1-h_width1, 0, 0, lc2};
Point(8) =  {posx3+posx1+h_width1, 0, 0, lc2};
Point(9) =  {posx3+posx1+h_width1, -d, 0, lc2};
Point(10) = {posx3+posx1-h_width1, -d, 0, lc2};
// Center
Point(15) =  {posx3+posx1, 0, 0, lc3};
Point(16) =  {posx3+posx1, -d, 0, lc3};

Point(11) = {posx2-h_width2, 0, 0, lc2};
Point(12) = {posx2+h_width2, 0, 0, lc2};
Point(13) = {posx2+h_width2, -d, 0, lc2};
Point(14) = {posx2-h_width2, -d, 0, lc2};
// Center
Point(17) = {posx2, 0, 0, lc4};
Point(18) = {posx2, -d, 0, lc4};



Line(1) = {1, 11};
Line(2) = {11, 17};
Line(3) = {17, 12};
Line(4) = {12, 5};
Line(5) = {5, 7};
Line(6) = {7, 15};
Line(7) = {15, 8};
Line(8) = {8, 4};
Line(9) = {4, 3};
Line(10) = {3, 9};
Line(11) = {9, 16};
Line(12) = {16, 10};
Line(13) = {10, 6};
Line(14) = {6, 13};
Line(15) = {13, 18};
Line(16) = {18, 14};
Line(17) = {14, 2};
Line(18) = {2, 1};
Line(19) = {11, 14};
Line(20) = {18, 17};
Line(21) = {12, 13};
Line(22) = {5, 6};
Line(23) = {7, 10};
Line(24) = {16, 15};
Line(25) = {8, 9};


Line Loop(26) = {1, 19, 17, 18};
Plane Surface(27) = {26};
Line Loop(28) = {20, -2, 19, -16};
Plane Surface(29) = {28};
Line Loop(30) = {15, 20, 3, 21};
Plane Surface(31) = {30};
Line Loop(32) = {4, 22, 14, -21};
Plane Surface(33) = {32};
Line Loop(34) = {13, -22, 5, 23};
Plane Surface(35) = {34};
Line Loop(36) = {6, -24, 12, -23};
Plane Surface(37) = {36};
Line Loop(38) = {11, 24, 7, 25};
Plane Surface(39) = {38};
Line Loop(40) = {8, 9, 10, -25};
Plane Surface(41) = {40};


Physical Line(42) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Line(43) = {9};
Physical Line(44) = {10, 11, 12, 13, 14, 15, 16, 17};
Physical Line(45) = {18};


Physical Surface(1) = {27, 33, 35, 41};
Physical Surface(2) = {29, 31};
Physical Surface(3) = {37, 39};

// Template mesh geometry file for L shape.

d = 1; // grating period
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;

// input parameters in nm
L_nm = 0;
W_nm = 0;
r = 0;

// normalized input parameters
L = L_nm/d_in_nm;
W = W_nm/d_in_nm;

// mesh finesse parameters
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // mesh at circular arcs

// unitary cell
Point(1) = {-0.5*d,-0.5*dy, 0, lc};
Point(2) = {0.5*d,-0.5*dy, 0, lc};
Point(3) = {0.5*d,0.5*dy, 0, lc};
Point(4) = {-0.5*d,0.5*dy, 0, lc};

// L points
Point(5) = {0.5*W*r-0.5*L,0-0.5*L, 0, lc2};
Point(6) = {L-0.5*W-0.5*L,0-0.5*L, 0, lc};
Point(7) = {L-0.5*W-0.5*L,0.5*W-0.5*L, 0, lc};
Point(8) = {L-0.5*W-0.5*L+0.5*Sqrt(2)*W/2,0.5*W-0.5*L - 0.5*Sqrt(2)*W/2, 0, lc2}; // lower right corner
Point(9) = {L-0.5*L,0.5*W-0.5*L, 0, lc2};
Point(10) = {L-0.5*W-0.5*L+0.5*Sqrt(2)*W/2,0.5*W-0.5*L + 0.5*Sqrt(2)*W/2, 0, lc2};
Point(11) = {L-0.5*W-0.5*L,W-0.5*L, 0, lc};
Point(12) = {W+0.5*W*r-0.5*L,W-0.5*L, 0, lc2};
Point(13) = {W+0.5*W*r-0.5*L,W+0.5*W*r-0.5*L, 0, lc};
Point(14) = {W-0.5*L,W+0.5*W*r-0.5*L, 0, lc2};
Point(15) = {W-0.5*L,L-0.5*W-0.5*L, 0, lc};
Point(16) = {0.5*W-0.5*L,L-0.5*W-0.5*L, 0, lc};
Point(17) = {0.5*W-0.5*L+0.5*Sqrt(2)*W/2,L-0.5*W-0.5*L+0.5*Sqrt(2)*W/2, 0, lc2}; // upper left corner
Point(18) = {0.5*W-0.5*L,L-0.5*L, 0, lc2};
Point(19) = {0.5*W-0.5*L-0.5*Sqrt(2)*W/2,L-0.5*W-0.5*L+0.5*Sqrt(2)*W/2, 0, lc2};
Point(20) = {0.0-0.5*L,L-0.5*W-0.5*L, 0, lc};
Point(21) = {0.0-0.5*L,0.5*W*r-0.5*L, 0, lc2};
Point(22) = {0.5*W*r-0.5*L,0.5*W*r-0.5*L, 0, lc};

// L lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Circle(6) = {6,7,8};
Circle(7) = {8,7,9};
Circle(8) = {9,7,10};
Circle(9) = {10,7,11};
Line(10) = {11,12};
Circle(11) = {12,13,14};
Line(12) = {14,15};
Circle(13) = {15,16,17};
Circle(14) = {17,16,18};
Circle(15) = {18,16,19};
Circle(16) = {19,16,20};
Line(17) = {20,21};
Circle(18) = {21,22,5};

// Line loops defining the surfaces
Line Loop(1) = {5,6,7,8,9,10,11,12,13,14,15,16,17,18};
Line Loop(2) = {1,2,3,4};

// Surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2,1};

// Physical Surfaces where materials get assigned
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Physical Surface(1) = {2}; //internal
Physical Surface(2) = {1};

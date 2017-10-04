// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;

// input geo params
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0; // horizontal rectangle side
a2 = 0; // vertical rectangle side
smooth = 0;

// input mesh params
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on edges
lc3 = lc/1; // on edges

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;

// derived params
r1 = (a1/(2*d_in_nm))*d;
r2 = (a2/(2*d_in_nm))*d;
s2 = r2*smooth;
r1r=r1-s2;
r2r=r2-s2;

// unit cell
Point(1) = {0, 0, 0, lc}; // first side
Point(2) = {0, -0.5*dy, 0, lc};
Point(3) = {0, -dy, 0, lc};
Point(4) = {0.5*d, -dy, 0,lc};
Point(5) = {d, -dy, 0,lc};
Point(6) = {d, -0.5*dy, 0,lc};
Point(7) = {d, 0, 0,lc};
Point(8) = {0.5*d, 0, 0,lc};
Point(9) = {0.5*d, -0.5*dy, 0,lc2};

// rectangle
Point(10) = {0.5*d+r1, -0.5*dy, 0,lc2}; // upper right
Point(11) = {0.5*d+r1, -0.5*dy+r2r, 0,lc3};
Point(12) = {0.5*d+r1r, -0.5*dy+r2r, 0,lc2};
Point(13) = {0.5*d+r1r, -0.5*dy+r2, 0,lc3};
Point(14) = {0.5*d, -0.5*dy+r2, 0,lc2};
Point(15) = {0.5*d-r1r, -0.5*dy+r2, 0,lc3}; // upper left
Point(16) = {0.5*d-r1r, -0.5*dy+r2r, 0,lc2};
Point(17) = {0.5*d-r1, -0.5*dy+r2r, 0,lc3};
Point(18) = {0.5*d-r1, -0.5*dy, 0,lc2};
Point(19) = {0.5*d-r1, -0.5*dy-r2r, 0,lc3}; //lower left
Point(20) = {0.5*d-r1r, -0.5*dy-r2r, 0,lc2};
Point(21) = {0.5*d-r1r, -0.5*dy-r2, 0,lc3};
Point(22) = {0.5*d, -0.5*dy-r2, 0,lc2};
Point(23) = {0.5*d+r1r, -0.5*dy-r2, 0,lc3}; //lower rigth
Point(24) = {0.5*d+r1r, -0.5*dy-r2r, 0,lc2};
Point(25) = {0.5*d+r1, -0.5*dy-r2r, 0,lc3};


// LINES: unit cell plus crossings
Line(1) = {1,2}; // perimeter
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {2,18}; // horizontal diameter
Line(10) = {18,9};
Line(11) = {9,10};
Line(12) = {10,6};
Line(13) = {8,14}; // vertical diameter
Line(14) = {14,9};
Line(15) = {9,22};
Line(16) = {22,4};

// LINES: rectangle
Line(21) = {10,11};
Ellipse(22) = {11,12,13,13};
Line(23) = {13,14};
Line(24) = {14,15};
Ellipse(25) = {15,16,17,17};
Line(26) = {17,18};
Line(27) = {18,19};
Ellipse(28) = {19,20,21,21};
Line(29) = {21,22};
Line(30) = {22,23};
Ellipse(31) = {23,24,25,25};
Line(32) = {25,10};

// LINELOOPS: external
Line Loop(41) = {6,7,13,-23,-22,-21,12};
Plane Surface(41) = {41};
Line Loop(42) = {-13,8,1,9,-26,-25,-24};
Plane Surface(42) = {42};
Line Loop(43) = {-16,-29,-28,-27,-9,2,3};
Plane Surface(43) = {43};
Line Loop(44) = {5,-12,-32,-31,-30,16,4};
Plane Surface(44) = {44};

// LINELOOPS: internal
Line Loop(45) = {32,-11,15,30,31};
Plane Surface(45) = {45};
Line Loop(46) = {21,22,23,14,11};
Plane Surface(46) = {46};
Line Loop(47) = {-14,24,25,26,10};
Plane Surface(47) = {47};
Line Loop(48) = {-15,-10,27,28,29};
Plane Surface(48) = {48};

//PHYSICAL ENTITIES
Physical Line(1) = {1,2}; // external
Physical Line(2) = {3,4};
Physical Line(3) = {5,6};
Physical Line(4) = {7,8};

Physical Surface(1) = {41,42,43,44}; //internal
Physical Surface(2) = {45,46,47,48};

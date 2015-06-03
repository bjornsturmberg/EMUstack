// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;
d_in_nm = 0; // input geo params
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0;
a2 = 0;
gap = 0;
ellipticity = 0;
lc = 0; // input mesh params
lc2 = lc/1; // on cylinder surfaces and center
lc3 = lc/1; // on gap
center1 = d_in_nm/2 - a1/2 - gap/2; // derived params
center2 = d_in_nm/2 + a2/2 + gap/2;
r1 = (a1/(2*d_in_nm))*d;
r2 = (a2/(2*d_in_nm))*d;
c1 = (center1/d_in_nm)*d;
c2 = (center2/d_in_nm)*d;
ell = ellipticity;

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;


// unit cell
Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// first circle
Point(5) = {-hx+c1, -hy/2, 0,lc2};
Point(6) = {-hx+c1, -hy/2+(r1-ell*r1), 0, lc2};
Point(7) = {-hx+c1-r1, -hy/2, 0, lc2};
Point(8) = {-hx+c1, -hy/2-(r1-ell*r1), 0, lc2};
Point(9) = {-hx+c1+r1, -hy/2, 0, lc3};

// second circle
Point(10) = {-hx+c2, -hy/2, 0,lc2};
Point(11) = {-hx+c2, -hy/2+(r2-ell*r2), 0, lc2};
Point(12) = {-hx+c2-r2, -hy/2, 0, lc3};
Point(13) = {-hx+c2, -hy/2-(r2-ell*r2), 0, lc2};
Point(14) = {-hx+c2+r2, -hy/2, 0, lc2};

// middle axis points
Point(15) = {-hx, -hy/2, 0, lc};
Point(16) = {-hx+d, -hy/2, 0, lc};

// now connecting the dots :)
Line(1) = {1,4}; // unit cell perimeter
Line(2) = {2,3};
Line(3) = {1,15};
Line(4) = {15,2};
Line(5) = {4,16};
Line(6) = {16,3};
Line(7) = {15,7}; // horizontal middle axis
Line(8) = {7,5};
Line(9) = {5,9};
Line(10) = {9,12};
Line(11) = {12,10};
Line(12) = {10,14};
Line(13) = {14,16};
Ellipsis(14) = {9,5,6,6}; // left circle
Ellipsis(15) = {6,5,7,7};
Ellipsis(16) = {7,5,8,8};
Ellipsis(17) = {8,5,9,9};
Ellipsis(18) = {14,10,11,11}; // right circle
Ellipsis(19) = {11,10,12,12};
Ellipsis(20) = {12,10,13,13};
Ellipsis(21) = {13,10,14,14};

// now the line loops to define the surfaces
// starting from the top right corner
Line Loop(22) = {5,-13,18,19,-10,14,15,-7,-3,1};
Plane Surface(23) = {22};
Line Loop(24) = {6,-2,-4,7,16,17,10,20,21,13};
Plane Surface(25) = {24};
Line Loop(26) = {-12,-11,-19,-18};
Plane Surface(27) = {26};
Line Loop(28) = {-9,-8,-15,-14};
Plane Surface(29) = {28};
Line Loop(30) = {-17,-16,8,9};
Plane Surface(31) = {30};
Line Loop(32) = {-21,-20,11,12};
Plane Surface(33) = {32};

// physical lines: what are they?
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3,4};
Physical Line(4) = {5,6};

// the real stuff
Physical Surface(1) = {23,25};
Physical Surface(2) = {27,29,31,33};

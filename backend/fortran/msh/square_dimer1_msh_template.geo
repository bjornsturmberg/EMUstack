// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;

// input geo params
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0;
a2 = 0;
gap = 0;
ellipticity = 0;
smooth = 0;

// input mesh params
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces and center
lc3 = lc/1; // on gap

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;

// derived params
center1 = d_in_nm/2 - a1/2 - gap/2;
center2 = d_in_nm/2 + a2/2 + gap/2;
r1 = (a1/(2*d_in_nm))*d;
r2 = (a2/(2*d_in_nm))*d;
c1 = (center1/d_in_nm)*d;
c2 = (center2/d_in_nm)*d;
ell = ellipticity;
s1 = r1*smooth;
s2 = r2*smooth;
r1r=r1-s1;
r2r=r2-s2;

// unit cell
Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// first square
Point(5) = {-hx+c1, -hy/2, 0,lc2}; // center
Point(6) = {-hx+c1+r1, -hy/2, 0, lc3}; // center right
Point(7) = {-hx+c1+r1r+s1, -hy/2+r1r, 0,lc3}; // top right corner
Point(8) = {-hx+c1+r1r, -hy/2+r1r, 0,lc2};
Point(9) = {-hx+c1+r1r, -hy/2+r1r+s1, 0,lc2};
Point(10) = {-hx+c1-r1r, -hy/2+r1r+s1, 0,lc2}; // top left corner
Point(11) = {-hx+c1-r1r, -hy/2+r1r, 0,lc2};
Point(12) = {-hx+c1-r1r-s1, -hy/2+r1r, 0,lc2};
Point(13) = {-hx+c1-r1, -hy/2, 0, lc2}; // center right
Point(14) = {-hx+c1-r1r-s1, -hy/2-r1r, 0,lc2}; // bottom left corner
Point(15) = {-hx+c1-r1r, -hy/2-r1r, 0,lc2};
Point(16) = {-hx+c1-r1r, -hy/2-r1r-s1, 0,lc2};
Point(17) = {-hx+c1+r1r, -hy/2-r1r-s1, 0,lc2}; // bottom right corner
Point(18) = {-hx+c1+r1r, -hy/2-r1r, 0,lc2};
Point(19) = {-hx+c1+r1r+s1, -hy/2-r1r, 0,lc3};

// second square
Point(20) = {-hx+c2, -hy/2, 0,lc2}; // center
Point(21) = {-hx+c2+r2, -hy/2, 0, lc2}; // center right
Point(22) = {-hx+c2+r2r+s2, -hy/2+r2r, 0,lc2}; // top right corner
Point(23) = {-hx+c2+r2r, -hy/2+r2r, 0,lc2};
Point(24) = {-hx+c2+r2r, -hy/2+r2r+s2, 0,lc2};
Point(25) = {-hx+c2-r2r, -hy/2+r2r+s2, 0,lc2}; // top left corner
Point(26) = {-hx+c2-r2r, -hy/2+r2r, 0,lc2};
Point(27) = {-hx+c2-r2r-s2, -hy/2+r2r, 0,lc3};
Point(28) = {-hx+c2-r2, -hy/2, 0, lc3}; // center right
Point(29) = {-hx+c2-r2r-s2, -hy/2-r2r, 0,lc3}; // bottom left corner
Point(30) = {-hx+c2-r2r, -hy/2-r2r, 0,lc2};
Point(31) = {-hx+c2-r2r, -hy/2-r2r-s2, 0,lc2};
Point(32) = {-hx+c2+r2r, -hy/2-r2r-s2, 0,lc2}; // bottom right corner
Point(33) = {-hx+c2+r2r, -hy/2-r2r, 0,lc2};
Point(34) = {-hx+c2+r2r+s2, -hy/2-r2r, 0,lc2};

// middle axis points
Point(35) = {-hx, -hy/2, 0,lc};
Point(36) = {-hx+d, -hy/2, 0,lc};

// now connecting the dots :)
Line(1) = {1,4}; // unit cell perimeter
Line(2) = {2,3};
Line(3) = {1,35};
Line(4) = {35,2};
Line(5) = {4,36};
Line(6) = {36,3};
Line(7) = {35,13}; // horizontal middle axis
Line(8) = {13,5};
Line(9) = {5,6};
Line(10) = {6,28};
Line(11) = {28,20};
Line(12) = {20,21};
Line(13) = {21,36};
Line(14) = {6,7}; // left square
Ellipsis(15) = {7,8,9,9};
Line(16) = {9,10};
Ellipsis(17) = {10,11,12,12};
Line(18) = {12,13};
Line(19) = {13,14};
Ellipsis(20) = {14,15,16,16};
Line(21) = {16,17};
Ellipsis(22) = {17,18,19,19};
Line(23) = {19,6};
Line(24) = {21,22}; // right square
Ellipsis(25) = {22,23,24,24};
Line(26) = {24,25};
Ellipsis(27) = {25,26,27,27};
Line(28) = {27,28};
Line(29) = {28,29};
Ellipsis(30) = {29,30,31,31};
Line(31) = {31,32};
Ellipsis(32) = {32,33,34,34};
Line(33) = {34,21};

// now the line loops to define the surfaces
// starting from the top right corner
Line Loop(34) = {5,-13,24,25,26,27,28,-10,14,15,16,17,18,-7,-3,1};
Plane Surface(35) = {34};
Line Loop(36) = {6,-2,-4,7,19,20,21,22,23,10,29,30,31,32,33,13};
Plane Surface(37) = {36};
Line Loop(38) = {-12,-11,-28,-27,-26,-25,-24};
Plane Surface(39) = {38};
Line Loop(40) = {-9,-8,-18,-17,-16,-15,-14};
Plane Surface(41) = {40};
Line Loop(42) = {-23,-22,-21,-20,-19,8,9};
Plane Surface(43) = {42};
Line Loop(44) = {-33,-32,-31,-30,-29,11,12};
Plane Surface(45) = {44};

// physical lines: what are they?
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3,4};
Physical Line(4) = {5,6};

// the real stuff
Physical Surface(1) = {35,37};
Physical Surface(2) = {39,41,43,45};

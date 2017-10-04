// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;

// input geo params
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0; // rect1 horizontal side
b1 = 0; // rect1 vertical side
t = 0;
ellipticity = 0;
smooth = 0;

// input mesh params
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces and center
lc3 = lc/1; // on gap

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;

// derived params
center1 = d_in_nm/2; // rectangle centers in nm

r1 = (a1/(2*d_in_nm))*d; // normalized inner rectangles horizontal sides
r1Big = ((a1+2*t)/(2*d_in_nm))*d; // normalized outer rectangles horizontal sides

r1b = (b1/(2*d_in_nm))*d; // normalized inner rectangles vertical sides
r1bBig = ((b1+2*t)/(2*d_in_nm))*d; // normalized outer rectangles vertical sides

c1 = (center1/d_in_nm)*d; // normalized centers
ell = ellipticity;

s1 = r1b*smooth; // normalized inner rectangles curvature radii
s1Big = r1bBig*smooth; // normalized outer rectangles curvature radii

r1r=r1-s1; // normalized inner rectangles horizontal inner points
r1rBig=r1Big-s1Big; // normalized outer rectangles horizontal inner points

r1rb=r1b-s1; // normalized inner rectangles vertical inner points
r1rbBig=r1bBig-s1Big; // normalized outer rectangles vertical inner points

// unit cell
Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// first square (inner)
Point(5) = {-hx+c1, -hy/2, 0,lc2}; // center
Point(6) = {-hx+c1+r1, -hy/2, 0, lc2}; // center right
Point(7) = {-hx+c1+r1r+s1, -hy/2+r1rb, 0,lc2}; // top right corner
Point(8) = {-hx+c1+r1r, -hy/2+r1rb, 0,lc2};
Point(9) = {-hx+c1+r1r, -hy/2+r1rb+s1, 0,lc2};
Point(10) = {-hx+c1-r1r, -hy/2+r1rb+s1, 0,lc2}; // top left corner
Point(11) = {-hx+c1-r1r, -hy/2+r1rb, 0,lc2};
Point(12) = {-hx+c1-r1r-s1, -hy/2+r1rb, 0,lc2};
Point(13) = {-hx+c1-r1, -hy/2, 0, lc2}; // center right
Point(14) = {-hx+c1-r1r-s1, -hy/2-r1rb, 0,lc2}; // bottom left corner
Point(15) = {-hx+c1-r1r, -hy/2-r1rb, 0,lc2};
Point(16) = {-hx+c1-r1r, -hy/2-r1rb-s1, 0,lc2};
Point(17) = {-hx+c1+r1r, -hy/2-r1rb-s1, 0,lc2}; // bottom right corner
Point(18) = {-hx+c1+r1r, -hy/2-r1rb, 0,lc2};
Point(19) = {-hx+c1+r1r+s1, -hy/2-r1rb, 0,lc2};


// first square (outer)
Point(106) = {-hx+c1+r1Big, -hy/2, 0, lc3}; // center right
Point(107) = {-hx+c1+r1rBig+s1Big, -hy/2+r1rbBig, 0,lc2}; // top right corner
Point(108) = {-hx+c1+r1rBig, -hy/2+r1rbBig, 0,lc2};
Point(109) = {-hx+c1+r1rBig, -hy/2+r1rbBig+s1Big, 0,lc2};
Point(110) = {-hx+c1-r1rBig, -hy/2+r1rbBig+s1Big, 0,lc2}; // top left corner
Point(111) = {-hx+c1-r1rBig, -hy/2+r1rbBig, 0,lc2};
Point(112) = {-hx+c1-r1rBig-s1Big, -hy/2+r1rbBig, 0,lc2};
Point(113) = {-hx+c1-r1Big, -hy/2, 0, lc2}; // center right
Point(114) = {-hx+c1-r1rBig-s1Big, -hy/2-r1rbBig, 0,lc2}; // bottom left corner
Point(115) = {-hx+c1-r1rBig, -hy/2-r1rbBig, 0,lc2};
Point(116) = {-hx+c1-r1rBig, -hy/2-r1rbBig-s1Big, 0,lc2};
Point(117) = {-hx+c1+r1rBig, -hy/2-r1rbBig-s1Big, 0,lc2}; // bottom right corner
Point(118) = {-hx+c1+r1rBig, -hy/2-r1rbBig, 0,lc2};
Point(119) = {-hx+c1+r1rBig+s1Big, -hy/2-r1rbBig, 0,lc2};

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
Line(7) = {35,113}; // horizontal middle axis
Line(8) = {113,13};
Line(9) = {13,5};
Line(10) = {5,6};
Line(11) = {6,106};
Line(12) = {106,36};

Line(114) = {6,7}; // left inner square
Ellipsis(115) = {7,8,9,9};
Line(116) = {9,10};
Ellipsis(117) = {10,11,12,12};
Line(118) = {12,13};
Line(119) = {13,14};
Ellipsis(120) = {14,15,16,16};
Line(121) = {16,17};
Ellipsis(122) = {17,18,19,19};
Line(123) = {19,6};

Line(214) = {106,107}; // left outer square
Ellipsis(215) = {107,108,109,109};
Line(216) = {109,110};
Ellipsis(217) = {110,111,112,112};
Line(218) = {112,113};
Line(219) = {113,114};
Ellipsis(220) = {114,115,116,116};
Line(221) = {116,117};
Ellipsis(222) = {117,118,119,119};
Line(223) = {119,106};

// now the line loops to define the surfaces
// starting from the top right corner
Line Loop(34) = {5,-12,214,215,216,217,218,-7,-3,1}; // upper big surface
Plane Surface(35) = {34};
Line Loop(36) = {6,-2,-4,7,219,220,221,222,223,12}; // lower big surface
Plane Surface(37) = {36};

Line Loop(38) = {114,115,116,117,118,9,10}; // upper left inner
Plane Surface(39) = {38};
Line Loop(40) = {-123,-122,-121,-120,-119,9,10}; // lower left inner
Plane Surface(41) = {40};
Line Loop(42) = {214,215,216,217,218,8,-118,-117,-116,-115,-114,11}; // upper left outer
Plane Surface(43) = {42};
Line Loop(44) = {-223,-222,-221,-220,-219,8,119,120,121,122,123,11}; // lower left outer
Plane Surface(45) = {44};

// physical lines: what are they?
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3,4};
Physical Line(4) = {5,6};

// the real stuff
Physical Surface(1) = {35,37}; // medium
Physical Surface(2) = {39,41}; // dimer
Physical Surface(3) = {43,45}; // dimer shell

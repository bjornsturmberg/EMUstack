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

// choosing long and short side
ss=a2;
ls=a1;
If(a1<a2)
ss=a1;
ls=a2;
EndIf

// derived params horizontal arm
r1 = (ls/(2*d_in_nm))*d;
r2 = (ss/(2*d_in_nm))*d;
s2 = r2*smooth;
r1r=r1-s2;
r2r=r2-s2;

// derived params vertical arm
r1v = (ss/(2*d_in_nm))*d;
r2v = (ls/(2*d_in_nm))*d;
s2v = r2*smooth;
r1rv=r1v-s2v;
r2rv=r2v-s2v;


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

// Center
Point(10) = {0.5*d+r2, -0.5*dy+r2, 0,lc2};
Point(11) = {0.5*d-r2, -0.5*dy+r2, 0,lc2};
Point(12) = {0.5*d-r2, -0.5*dy-r2, 0,lc2};
Point(13) = {0.5*d+r2, -0.5*dy-r2, 0,lc2};

// Inner mesh points
Point(50) = {0.5*d+0.8*r2, -0.5*dy+0.8*r2, 0,lc2};
Point(51) = {0.5*d-0.8*r2, -0.5*dy+0.8*r2, 0,lc2};
Point(52) = {0.5*d-0.8*r2, -0.5*dy-0.8*r2, 0,lc2};
Point(53) = {0.5*d+0.8*r2, -0.5*dy-0.8*r2, 0,lc2};

// horizontal rectangle
Point(110) = {0.5*d+r1, -0.5*dy, 0,lc2}; // upper right
Point(111) = {0.5*d+r1, -0.5*dy+r2r, 0,lc3};
Point(112) = {0.5*d+r1r, -0.5*dy+r2r, 0,lc2};
Point(113) = {0.5*d+r1r, -0.5*dy+r2, 0,lc3};
Point(114) = {0.5*d-r1r, -0.5*dy+r2, 0,lc3}; // upper left
Point(115) = {0.5*d-r1r, -0.5*dy+r2r, 0,lc2};
Point(116) = {0.5*d-r1, -0.5*dy+r2r, 0,lc3};
Point(117) = {0.5*d-r1, -0.5*dy, 0,lc2};
Point(118) = {0.5*d-r1, -0.5*dy-r2r, 0,lc3}; //lower left
Point(119) = {0.5*d-r1r, -0.5*dy-r2r, 0,lc2};
Point(120) = {0.5*d-r1r, -0.5*dy-r2, 0,lc3};
Point(121) = {0.5*d+r1r, -0.5*dy-r2, 0,lc3}; //lower rigth
Point(122) = {0.5*d+r1r, -0.5*dy-r2r, 0,lc2};
Point(123) = {0.5*d+r1, -0.5*dy-r2r, 0,lc3};

// vertical rectangle
Point(210) = {0.5*d+r1v, -0.5*dy+r2rv, 0,lc3};
Point(211) = {0.5*d+r1rv, -0.5*dy+r2rv, 0,lc2};
Point(212) = {0.5*d+r1rv, -0.5*dy+r2v, 0,lc3};
Point(213) = {0.5*d, -0.5*dy+r2v, 0,lc2};
Point(214) = {0.5*d-r1rv, -0.5*dy+r2v, 0,lc3}; // upper left
Point(215) = {0.5*d-r1rv, -0.5*dy+r2rv, 0,lc2};
Point(216) = {0.5*d-r1v, -0.5*dy+r2rv, 0,lc3};
Point(217) = {0.5*d-r1v, -0.5*dy-r2rv, 0,lc3}; //lower left
Point(218) = {0.5*d-r1rv, -0.5*dy-r2rv, 0,lc2};
Point(219) = {0.5*d-r1rv, -0.5*dy-r2v, 0,lc3};
Point(220) = {0.5*d, -0.5*dy-r2v, 0,lc2};
Point(221) = {0.5*d+r1rv, -0.5*dy-r2v, 0,lc3}; //lower rigth
Point(222) = {0.5*d+r1rv, -0.5*dy-r2rv, 0,lc2};
Point(223) = {0.5*d+r1v, -0.5*dy-r2rv, 0,lc3};

// LINES: unit cell plus crossings
Line(1) = {1,2}; // perimeter
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {2,117}; // horizontal diameter
Line(10) = {117,9};
Line(11) = {9,110};
Line(12) = {110,6};
Line(13) = {8,213}; // vertical diameter
Line(14) = {213,9};
Line(15) = {9,220};
Line(16) = {220,4};

// LINES: cross
Line(21) = {110,111};
Ellipse(22) = {111,112,113,113};
Line(23) = {113,10};
Line(24) = {10,210};
Ellipse(25) = {210,211,212,212};
Line(26) = {212,213};
Line(27) = {213,214};
Ellipse(28) = {214,215,216,216};
Line(29) = {216,11};
Line(30) = {11,114};
Ellipse(31) = {114,115,116,116};
Line(32) = {116,117};
Line(33) = {117,118};
Ellipse(34) = {118,119,120,120};
Line(35) = {120,12};
Line(36) = {12,217};
Ellipse(37) = {217,218,219,219};
Line(38) = {219,220};
Line(39) = {220,221};
Ellipse(40) = {221,222,223,223};
Line(41) = {223,13};
Line(42) = {13,121};
Ellipse(43) = {121,122,123,123};
Line(44) = {123,110};

// LINELOOPS: external
Line Loop(41) = {6,7,13,-26,-25,-24,-23,-22,-21,12};
Plane Surface(41) = {41};
Line Loop(42) = {-13,8,1,9,-32,-31,-30,-29,-28,-27};
Plane Surface(42) = {42};
Line Loop(43) = {-16,-38,-37,-36,-35,-34,-33,-9,2,3};
Plane Surface(43) = {43};
Line Loop(44) = {5,-12,-44,-43,-42,-41,-40,-39,16,4};
Plane Surface(44) = {44};

// LINELOOPS: internal
Line Loop(45) = {21,22,23,24,25,26,14,11};
Plane Surface(45) = {45};
Line Loop(46) = {27,28,29,30,31,32,10,-14};
Plane Surface(46) = {46};
Line Loop(47) = {33,34,35,36,37,38,-15,-10};
Plane Surface(47) = {47};
Line Loop(48) = {39,40,41,42,43,44,-11,15};
Plane Surface(48) = {48};

Point{50} In Surface{45};
Point{51} In Surface{46};
Point{52} In Surface{47};
Point{53} In Surface{48};

//PHYSICAL ENTITIES
Physical Line(1) = {1,2}; // external
Physical Line(2) = {3,4};
Physical Line(3) = {5,6};
Physical Line(4) = {7,8};

Physical Surface(1) = {41,42,43,44}; //internal
Physical Surface(2) = {45,46,47,48};

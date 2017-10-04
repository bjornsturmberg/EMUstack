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
a2 = 0; // rect2 horizontal side
b1 = 0; // rect1 vertical side
b2 = 0; // rect2 vertical side
gap = 0;
thickness = 0;
ellipticity = 0;
smooth = 0;

// input mesh params
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces and center
lc3 = lc/1; // on gap

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;

// derived params
center1 = d_in_nm/2 - a1/2 - gap/2; // rectangle centers in nm
center2 = d_in_nm/2 + a2/2 + gap/2;

r1 = (a1/(2*d_in_nm))*d; // normalized inner rectangles horizontal sides
r2 = (a2/(2*d_in_nm))*d;
r1Big = ((a1+2*t)/(2*d_in_nm))*d; // normalized outer rectangles horizontal sides
r2Big = ((a2+2*t)/(2*d_in_nm))*d;

r1b = (b1/(2*d_in_nm))*d; // normalized inner rectangles vertical sides
r2b = (b2/(2*d_in_nm))*d;
r1bBig = ((b1+2*t)/(2*d_in_nm))*d; // normalized outer rectangles vertical sides
r2bBig = ((b2+2*t)/(2*d_in_nm))*d;

c1 = (center1/d_in_nm)*d; // normalized centers
c2 = (center2/d_in_nm)*d;
ell = ellipticity;

s1 = r1b*smooth; // normalized inner rectangles curvature radii
s2 = r2b*smooth;
s1Big = r1bBig*smooth; // normalized outer rectangles curvature radii
s2Big = r2bBig*smooth;

r1r=r1-s1; // normalized inner rectangles horizontal inner points
r2r=r2-s2;
r1rBig=r1Big-s1Big; // normalized outer rectangles horizontal inner points
r2rBig=r2Big-s2Big;

r1rb=r1b-s1; // normalized inner rectangles vertical inner points
r2rb=r2b-s2;
r1rbBig=r1bBig-s1Big; // normalized outer rectangles vertical inner points
r2rbBig=r2bBig-s2Big;

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

// second square (inner)
Point(20) = {-hx+c2, -hy/2, 0,lc2}; // center
Point(21) = {-hx+c2+r2, -hy/2, 0, lc2}; // center right
Point(22) = {-hx+c2+r2r+s2, -hy/2+r2rb, 0,lc2}; // top right corner
Point(23) = {-hx+c2+r2r, -hy/2+r2rb, 0,lc2};
Point(24) = {-hx+c2+r2r, -hy/2+r2rb+s2, 0,lc2};
Point(25) = {-hx+c2-r2r, -hy/2+r2rb+s2, 0,lc2}; // top left corner
Point(26) = {-hx+c2-r2r, -hy/2+r2rb, 0,lc2};
Point(27) = {-hx+c2-r2r-s2, -hy/2+r2rb, 0,lc2};
Point(28) = {-hx+c2-r2, -hy/2, 0, lc3}; // center right
Point(29) = {-hx+c2-r2r-s2, -hy/2-r2rb, 0,lc2}; // bottom left corner
Point(30) = {-hx+c2-r2r, -hy/2-r2rb, 0,lc2};
Point(31) = {-hx+c2-r2r, -hy/2-r2rb-s2, 0,lc2};
Point(32) = {-hx+c2+r2r, -hy/2-r2rb-s2, 0,lc2}; // bottom right corner
Point(33) = {-hx+c2+r2r, -hy/2-r2rb, 0,lc2};
Point(34) = {-hx+c2+r2r+s2, -hy/2-r2rb, 0,lc2};

// second square (outer)
Point(121) = {-hx+c2+r2Big, -hy/2, 0, lc2}; // center right
Point(122) = {-hx+c2+r2rBig+s2Big, -hy/2+r2rbBig, 0,lc2}; // top right corner
Point(123) = {-hx+c2+r2rBig, -hy/2+r2rbBig, 0,lc2};
Point(124) = {-hx+c2+r2rBig, -hy/2+r2rbBig+s2Big, 0,lc2};
Point(125) = {-hx+c2-r2rBig, -hy/2+r2rbBig+s2Big, 0,lc2}; // top left corner
Point(126) = {-hx+c2-r2rBig, -hy/2+r2rbBig, 0,lc2};
Point(127) = {-hx+c2-r2rBig-s2Big, -hy/2+r2rbBig, 0,lc2};
Point(128) = {-hx+c2-r2Big, -hy/2, 0, lc3}; // center right
Point(129) = {-hx+c2-r2rBig-s2Big, -hy/2-r2rbBig, 0,lc2}; // bottom left corner
Point(130) = {-hx+c2-r2rBig, -hy/2-r2rbBig, 0,lc2};
Point(131) = {-hx+c2-r2rBig, -hy/2-r2rbBig-s2Big, 0,lc2};
Point(132) = {-hx+c2+r2rBig, -hy/2-r2rbBig-s2Big, 0,lc2}; // bottom right corner
Point(133) = {-hx+c2+r2rBig, -hy/2-r2rbBig, 0,lc2};
Point(134) = {-hx+c2+r2rBig+s2Big, -hy/2-r2rbBig, 0,lc2};

// middle axis points
Point(35) = {-hx, -hy/2, 0,lc};
Point(36) = {-hx+d, -hy/2, 0,lc};


// If gap larger than zero
If(gap-2.0*t>0)

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
Line(12) = {106,128};
Line(13) = {128,28};
Line(14) = {28,20};
Line(15) = {20,21};
Line(16) = {21,121};
Line(17) = {121,36};

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

Line(124) = {21,22}; // right inner square
Ellipsis(125) = {22,23,24,24};
Line(126) = {24,25};
Ellipsis(127) = {25,26,27,27};
Line(128) = {27,28};
Line(129) = {28,29};
Ellipsis(130) = {29,30,31,31};
Line(131) = {31,32};
Ellipsis(132) = {32,33,34,34};
Line(133) = {34,21};

Line(224) = {121,122}; // right outer square
Ellipsis(225) = {122,123,124,124};
Line(226) = {124,125};
Ellipsis(227) = {125,126,127,127};
Line(228) = {127,128};
Line(229) = {128,129};
Ellipsis(230) = {129,130,131,131};
Line(231) = {131,132};
Ellipsis(232) = {132,133,134,134};
Line(233) = {134,121};

// now the line loops to define the surfaces
// starting from the top right corner
Line Loop(34) = {5,-17,224,225,226,227,228,-12,214,215,216,217,218,-7,-3,1}; // upper big surface
Plane Surface(35) = {34};
Line Loop(36) = {6,-2,-4,7,219,220,221,222,223,12,229,230,231,232,233,17}; // lower big surface
Plane Surface(37) = {36};

Line Loop(38) = {114,115,116,117,118,9,10}; // upper left inner
Plane Surface(39) = {38};
Line Loop(40) = {-123,-122,-121,-120,-119,9,10}; // lower left inner
Plane Surface(41) = {40};
Line Loop(42) = {214,215,216,217,218,8,-118,-117,-116,-115,-114,11}; // upper left outer
Plane Surface(43) = {42};
Line Loop(44) = {-223,-222,-221,-220,-219,8,119,120,121,122,123,11}; // lower left outer
Plane Surface(45) = {44};

Line Loop(46) = {124,125,126,127,128,14,15}; // upper right inner
Plane Surface(47) = {46};
Line Loop(48) = {-133,-132,-131,-130,-129,14,15}; // lower right inner
Plane Surface(49) = {48};
Line Loop(50) = {224,225,226,227,228,13,-128,-127,-126,-125,-124,16}; // upper right outer
Plane Surface(51) = {50};
Line Loop(52) = {-233,-232,-231,-230,-229,13,129,130,131,132,133,16}; // lower right outer
Plane Surface(53) = {52};

// physical lines: what are they?
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3,4};
Physical Line(4) = {5,6};

// the real stuff
Physical Surface(1) = {35,37}; // medium
Physical Surface(2) = {39,41,47,49}; // dimer
Physical Surface(3) = {43,45,51,53}; // dimer shell

EndIf


// If gap equal or smaller than zero
If(gap-2.0*t<=0)

// Removing middle points
Delete {Point {106,107,108,109,117,118,119};}
Delete {Point {125,126,127,128,129,130,131};}

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
Line(11) = {6,28};
Line(12) = {28,20};
Line(13) = {20,21};
Line(14) = {21,121};
Line(15) = {121,36};

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

Line(124) = {21,22}; // right inner square
Ellipsis(125) = {22,23,24,24};
Line(126) = {24,25};
Ellipsis(127) = {25,26,27,27};
Line(128) = {27,28};
Line(129) = {28,29};
Ellipsis(130) = {29,30,31,31};
Line(131) = {31,32};
Ellipsis(132) = {32,33,34,34};
Line(133) = {34,21};

Line(214) = {121,122}; // left outer square
Ellipsis(215) = {122,123,124,124};
Line(216) = {124,110};
Ellipsis(217) = {110,111,112,112};
Line(218) = {112,113};
Line(219) = {113,114};
Ellipsis(220) = {114,115,116,116};
Line(221) = {116,132};
Ellipsis(222) = {132,133,134,134};
Line(223) = {134,121};


// now the line loops to define the surfaces
// starting from the top right corner
Line Loop(34) = {5,-15,214,215,216,217,218,-7,-3,1}; // upper big surface
Plane Surface(35) = {34};
Line Loop(36) = {6,-2,-4,7,219,220,221,222,223,15}; // lower big surface
Plane Surface(37) = {36};

Line Loop(38) = {114,115,116,117,118,9,10}; // upper left inner
Plane Surface(39) = {38};
Line Loop(40) = {-123,-122,-121,-120,-119,9,10}; // lower left inner
Plane Surface(41) = {40};

Line Loop(46) = {124,125,126,127,128,12,13}; // upper right inner
Plane Surface(47) = {46};
Line Loop(48) = {-133,-132,-131,-130,-129,12,13}; // lower right inner
Plane Surface(49) = {48};
Line Loop(50) = {214,215,216,217,218,8,-118,-117,-116,-115,-114,11,-128,-127,-126,-125,-124,14}; // upper outer
Plane Surface(51) = {50};
Line Loop(52) = {-223,-222,-221,-220,-219,8,119,120,121,122,123,11,129,130,131,132,133,14}; // lower outer
Plane Surface(53) = {52};

// physical lines: what are they?
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3,4};
Physical Line(4) = {5,6};

// the real stuff
Physical Surface(1) = {35,37}; // medium
Physical Surface(2) = {39,41,47,49}; // dimer
Physical Surface(3) = {51,53}; // dimer shell

EndIf

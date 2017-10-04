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
t = 0;

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

// derived params horizontal arm Big
r1Big = ((ls+2*t)/(2*d_in_nm))*d;
r2Big = ((ss+2*t)/(2*d_in_nm))*d;
s2Big = r2Big*smooth;
r1rBig=r1Big-s2Big;
r2rBig=r2Big-s2Big;

// derived params vertical arm
r1v = (ss/(2*d_in_nm))*d;
r2v = (ls/(2*d_in_nm))*d;
s2v = r2*smooth;
r1rv=r1v-s2v;
r2rv=r2v-s2v;

// derived params vertical arm Big
r1vBig = ((ss+2*t)/(2*d_in_nm))*d;
r2vBig = ((ls+2*t)/(2*d_in_nm))*d;
s2vBig = r2Big*smooth;
r1rvBig=r1vBig-s2vBig;
r2rvBig=r2vBig-s2vBig;


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

// Center Big
Point(1010) = {0.5*d+r2Big, -0.5*dy+r2Big, 0,lc2};
Point(1011) = {0.5*d-r2Big, -0.5*dy+r2Big, 0,lc2};
Point(1012) = {0.5*d-r2Big, -0.5*dy-r2Big, 0,lc2};
Point(1013) = {0.5*d+r2Big, -0.5*dy-r2Big, 0,lc2};

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

// horizontal rectangle Big
Point(1110) = {0.5*d+r1Big, -0.5*dy, 0,lc2}; // upper right
Point(1111) = {0.5*d+r1Big, -0.5*dy+r2rBig, 0,lc3};
Point(1112) = {0.5*d+r1rBig, -0.5*dy+r2rBig, 0,lc2};
Point(1113) = {0.5*d+r1rBig, -0.5*dy+r2Big, 0,lc3};
Point(1114) = {0.5*d-r1rBig, -0.5*dy+r2Big, 0,lc3}; // upper left
Point(1115) = {0.5*d-r1rBig, -0.5*dy+r2rBig, 0,lc2};
Point(1116) = {0.5*d-r1Big, -0.5*dy+r2rBig, 0,lc3};
Point(1117) = {0.5*d-r1Big, -0.5*dy, 0,lc2};
Point(1118) = {0.5*d-r1Big, -0.5*dy-r2rBig, 0,lc3}; //lower left
Point(1119) = {0.5*d-r1rBig, -0.5*dy-r2rBig, 0,lc2};
Point(1120) = {0.5*d-r1rBig, -0.5*dy-r2Big, 0,lc3};
Point(1121) = {0.5*d+r1rBig, -0.5*dy-r2Big, 0,lc3}; //lower rigth
Point(1122) = {0.5*d+r1rBig, -0.5*dy-r2rBig, 0,lc2};
Point(1123) = {0.5*d+r1Big, -0.5*dy-r2rBig, 0,lc3};

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

// vertical rectangle Big
Point(1210) = {0.5*d+r1vBig, -0.5*dy+r2rvBig, 0,lc3};
Point(1211) = {0.5*d+r1rvBig, -0.5*dy+r2rvBig, 0,lc2};
Point(1212) = {0.5*d+r1rvBig, -0.5*dy+r2vBig, 0,lc3};
Point(1213) = {0.5*d, -0.5*dy+r2vBig, 0,lc2};
Point(1214) = {0.5*d-r1rvBig, -0.5*dy+r2vBig, 0,lc3}; // upper left
Point(1215) = {0.5*d-r1rvBig, -0.5*dy+r2rvBig, 0,lc2};
Point(1216) = {0.5*d-r1vBig, -0.5*dy+r2rvBig, 0,lc3};
Point(1217) = {0.5*d-r1vBig, -0.5*dy-r2rvBig, 0,lc3}; //lower left
Point(1218) = {0.5*d-r1rvBig, -0.5*dy-r2rvBig, 0,lc2};
Point(1219) = {0.5*d-r1rvBig, -0.5*dy-r2vBig, 0,lc3};
Point(1220) = {0.5*d, -0.5*dy-r2vBig, 0,lc2};
Point(1221) = {0.5*d+r1rvBig, -0.5*dy-r2vBig, 0,lc3}; //lower rigth
Point(1222) = {0.5*d+r1rvBig, -0.5*dy-r2rvBig, 0,lc2};
Point(1223) = {0.5*d+r1vBig, -0.5*dy-r2rvBig, 0,lc3};

// LINES: unit cell plus crossings
Line(1) = {1,2}; // perimeter
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {2,1117}; // horizontal diameter
Line(10) = {1117,117};
Line(11) = {117,9};
Line(12) = {9,110};
Line(13) = {110,1110};
Line(14) = {1110,6};
Line(15) = {8,1213}; // vertical diameter
Line(16) = {1213,213};
Line(17) = {213,9};
Line(18) = {9,220};
Line(19) = {220,1220};
Line(20) = {1220,4};

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

// LINES: cross Big
Line(1021) = {1110,1111};
Ellipse(1022) = {1111,1112,1113,1113};
Line(1023) = {1113,1010};
Line(1024) = {1010,1210};
Ellipse(1025) = {1210,1211,1212,1212};
Line(1026) = {1212,1213};
Line(1027) = {1213,1214};
Ellipse(1028) = {1214,1215,1216,1216};
Line(1029) = {1216,1011};
Line(1030) = {1011,1114};
Ellipse(1031) = {1114,1115,1116,1116};
Line(1032) = {1116,1117};
Line(1033) = {1117,1118};
Ellipse(1034) = {1118,1119,1120,1120};
Line(1035) = {1120,1012};
Line(1036) = {1012,1217};
Ellipse(1037) = {1217,1218,1219,1219};
Line(1038) = {1219,1220};
Line(1039) = {1220,1221};
Ellipse(1040) = {1221,1222,1223,1223};
Line(1041) = {1223,1013};
Line(1042) = {1013,1121};
Ellipse(1043) = {1121,1122,1123,1123};
Line(1044) = {1123,1110};

// LINELOOPS: external
Line Loop(41) = {6,7,15,-1026,-1025,-1024,-1023,-1022,-1021,14};
Plane Surface(41) = {41};
Line Loop(42) = {-15,8,1,9,-1032,-1031,-1030,-1029,-1028,-1027};
Plane Surface(42) = {42};
Line Loop(43) = {-20,-1038,-1037,-1036,-1035,-1034,-1033,-9,2,3};
Plane Surface(43) = {43};
Line Loop(44) = {5,-14,-1044,-1043,-1042,-1041,-1040,-1039,20,4};
Plane Surface(44) = {44};

// LINELOOPS: shell
Line Loop(45) = {1021,1022,1023,1024,1025,1026,16,-26,-25,-24,-23,-22,-21,13};
Plane Surface(45) = {45};
Line Loop(46) = {-16,1027,1028,1029,1030,1031,1032,10,-32,-31,-30,-29,-28,-27};
Plane Surface(46) = {46};
Line Loop(47) = {1033,1034,1035,1036,1037,1038,-19,-38,-37,-36,-35,-34,-33,-10};
Plane Surface(47) = {47};
Line Loop(48) = {19,1039,1040,1041,1042,1043,1044,-13,-44,-43,-42,-41,-40,-39};
Plane Surface(48) = {48};

// LINELOOPS: core
Line Loop(49) = {21,22,23,24,25,26,17,12};
Plane Surface(49) = {49};
Line Loop(50) = {27,28,29,30,31,32,11,-17};
Plane Surface(50) = {50};
Line Loop(51) = {33,34,35,36,37,38,-18,-11};
Plane Surface(51) = {51};
Line Loop(52) = {44,-12,18,39,40,41,42,43};
Plane Surface(52) = {52};

Point{50} In Surface{49};
Point{51} In Surface{50};
Point{52} In Surface{51};
Point{53} In Surface{52};

//PHYSICAL ENTITIES
Physical Line(1) = {1,2}; // external
Physical Line(2) = {3,4};
Physical Line(3) = {5,6};
Physical Line(4) = {7,8};

// the real stuff
Physical Surface(1) = {41,42,43,44}; // medium
Physical Surface(2) = {49,50,51,52}; // dimer
Physical Surface(3) = {45,46,47,48}; // dimer shell

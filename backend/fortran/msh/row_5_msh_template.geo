// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

d = 1; // grating period
ff = 0;
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 10;
radius1 = (a1/(2*d_in_nm))*d;
sep = 5;
sep_norm = sep/d_in_nm + radius1;
ellipticity = 0;
square = 0;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres
lc4 = lc/1; // centres of top and bottom

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Vertices of the triangles

Point(5) = {-hx+d/2, -hy/2, 0,lc3};
Point(6) = {-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(7) = {-hx+d/2-radius1, -hy/2, 0, lc2};
Point(8) = {-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(9) = {-hx+d/2+radius1, -hy/2, 0, lc2};
Point(10) = {-hx+d/2, 0, 0, lc4};
Point(11) = {0,-hy/2, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc4};
Point(13) = {d, -hy/2, 0, lc};

Point(14) = {radius1+sep_norm-hx+d/2, -hy/2, 0,lc3};
Point(15) = {radius1+sep_norm-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(16) = {radius1+sep_norm-hx+d/2-radius1, -hy/2, 0, lc2};
Point(17) = {radius1+sep_norm-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(18) = {radius1+sep_norm-hx+d/2+radius1, -hy/2, 0, lc2};
Point(19) = {radius1+sep_norm-hx+d/2, 0, 0, lc4};
Point(20) = {radius1+sep_norm-hx+d/2, -hy, 0, lc4};

Point(21) = {-radius1-sep_norm-hx+d/2, -hy/2, 0,lc3};
Point(22) = {-radius1-sep_norm-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(23) = {-radius1-sep_norm-hx+d/2-radius1, -hy/2, 0, lc2};
Point(24) = {-radius1-sep_norm-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(25) = {-radius1-sep_norm-hx+d/2+radius1, -hy/2, 0, lc2};
Point(26) = {-radius1-sep_norm-hx+d/2, 0, 0, lc4};
Point(27) = {-radius1-sep_norm-hx+d/2, -hy, 0, lc4};

Point(28) = {2*radius1+2*sep_norm-hx+d/2, -hy/2, 0,lc3};
Point(29) = {2*radius1+2*sep_norm-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(30) = {2*radius1+2*sep_norm-hx+d/2-radius1, -hy/2, 0, lc2};
Point(31) = {2*radius1+2*sep_norm-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(32) = {2*radius1+2*sep_norm-hx+d/2+radius1, -hy/2, 0, lc2};
Point(33) = {2*radius1+2*sep_norm-hx+d/2, 0, 0, lc4};
Point(34) = {2*radius1+2*sep_norm-hx+d/2, -hy, 0, lc4};

Point(35) = {-2*radius1-2*sep_norm-hx+d/2, -hy/2, 0,lc3};
Point(36) = {-2*radius1-2*sep_norm-hx+d/2, -hy/2+(radius1-ellipticity*radius1), 0, lc2};
Point(37) = {-2*radius1-2*sep_norm-hx+d/2-radius1, -hy/2, 0, lc2};
Point(38) = {-2*radius1-2*sep_norm-hx+d/2, -hy/2-(radius1-ellipticity*radius1), 0, lc2};
Point(39) = {-2*radius1-2*sep_norm-hx+d/2+radius1, -hy/2, 0, lc2};
Point(40) = {-2*radius1-2*sep_norm-hx+d/2, 0, 0, lc4};
Point(41) = {-2*radius1-2*sep_norm-hx+d/2, -hy, 0, lc4};


Line(1) = {1,40};
Line(2) = {10,19};
Line(3) = {2,41};
Line(4) = {12,20};
Line(5) = {1,11};
Line(6) = {11,2};
Line(7) = {4,13};
Line(8) = {13,3};
Line(10) = {7,5};
Line(11) = {5,9};
Line(14) = {6,5};
Line(15) = {5,8};
Line(57) = {19, 33};
Line(58) = {34, 20};

Line(21) = {11, 37};
Line(22) = {10, 6};
Line(23) = {8, 12};
Line(24) = {9, 16};
Line(25) = {18, 30};
Line(26) = {19, 15};
Line(27) = {17, 20};
Line(32) = {15, 14};
Line(33) = {14, 17};
Line(34) = {14, 18};
Line(35) = {14, 16};

Physical Line(59) = {1, 2, 57};
Physical Line(60) = {7, 8};
Physical Line(61) = {58, 4, 3};
Physical Line(62) = {6, 5};

//If(ellipticity == 0)
//Circle(17) = {9,5,6};
//Circle(18) = {6,5,7};
//Circle(19) = {7,5,8};
//Circle(20) = {8,5,9};
//EndIf
If(square == 0)
Ellipsis(17) = {9,5,6,6};
Ellipsis(18) = {6,5,7,7};
Ellipsis(19) = {7,5,8,8};
Ellipsis(20) = {8,5,9,9};
Ellipse(28) = {18, 14, 15, 15};
Ellipse(29) = {18, 14, 17, 17};
Ellipse(30) = {17, 14, 16, 16};
Ellipse(31) = {16, 14, 15, 15};

Line(63) = {40, 26};
Line(64) = {26, 10};
Line(65) = {26, 22};
Line(66) = {22, 21};
Line(67) = {25, 7};
Line(68) = {25, 21};
Line(69) = {21, 24};
Line(70) = {27, 12};
Line(71) = {41, 27};
Line(72) = {35, 39};
Line(73) = {37, 35};
Line(74) = {23, 21};
Line(75) = {27, 24};
Line(76) = {38, 41};
Line(77) = {36, 35};
Line(78) = {38, 35};
Line(79) = {39, 23};
Line(80) = {36, 40};
Line(81) = {33, 4};
Line(82) = {34, 3};
Line(83) = {34, 31};
Line(84) = {29, 28};
Line(85) = {31, 28};
Line(86) = {33, 29};
Ellipse(87) = {37, 35, 36, 36};
Ellipse(88) = {39, 35, 36, 36};
Ellipse(89) = {39, 35, 38, 38};
Ellipse(90) = {38, 35, 37, 37};
Ellipse(91) = {23, 21, 22, 22};
Ellipse(92) = {25, 21, 22, 22};
Ellipse(93) = {25, 21, 24, 24};
Ellipse(94) = {24, 21, 23, 23};
Ellipse(95) = {30, 28, 29, 29};
Ellipse(96) = {29, 28, 32, 32};
Ellipse(97) = {32, 28, 31, 31};
Ellipse(98) = {31, 28, 30, 30};

Line Loop(99) = {5, 21, 87, 80, -1};
Plane Surface(100) = {99};
Line Loop(101) = {21, -90, 76, -3, -6};
Plane Surface(102) = {101};
Line Loop(103) = {76, 71, 75, 94, -79, 89};
Plane Surface(104) = {103};
Line Loop(105) = {75, -93, 67, 19, 23, -70};
Plane Surface(106) = {105};
Line Loop(107) = {23, 4, -27, 30, -24, -20};
Line Loop(108) = {79, 91, -65, -63, -80, -88};
Plane Surface(109) = {108};
Line Loop(110) = {65, -92, 67, -18, -22, -64};
Plane Surface(111) = {110};
Line Loop(112) = {22, -17, 24, 31, -26, -2};
Plane Surface(113) = {112};
Plane Surface(114) = {107};
Line(115) = {30, 28};
Line(116) = {28, 32};
Line(117) = {32, 13};
Line Loop(118) = {29, 27, -58, 83, 98, -25};
Plane Surface(119) = {118};
Line Loop(120) = {26, -28, 25, 95, -86, -57};
Plane Surface(121) = {120};
Line Loop(122) = {81, 7, -117, -96, -86};
Plane Surface(123) = {122};
Line Loop(124) = {97, -83, 82, -8, -117};
Plane Surface(125) = {124};
Line Loop(126) = {87, 77, -73};
Plane Surface(127) = {126};
Line Loop(128) = {73, -78, 90};
Plane Surface(129) = {128};
Line Loop(130) = {78, 72, 89};
Plane Surface(131) = {130};
Line Loop(132) = {72, 88, 77};
Plane Surface(133) = {132};
Line Loop(134) = {91, 66, -74};
Plane Surface(135) = {134};
Line Loop(136) = {74, 69, 94};
Plane Surface(137) = {136};
Line Loop(138) = {69, -93, 68};
Plane Surface(139) = {138};
Line Loop(140) = {68, -66, -92};
Plane Surface(141) = {140};
Line Loop(142) = {18, 10, -14};
Plane Surface(143) = {142};
Line Loop(144) = {10, 15, -19};
Plane Surface(145) = {144};
Line Loop(146) = {15, 20, -11};
Plane Surface(147) = {146};
Line Loop(148) = {11, 17, 14};
Plane Surface(149) = {148};
Line Loop(150) = {31, 32, 35};
Plane Surface(151) = {150};
Line Loop(152) = {35, -30, -33};
Plane Surface(153) = {152};
Line Loop(154) = {33, -29, -34};
Plane Surface(155) = {154};
Line Loop(156) = {34, 28, 32};
Plane Surface(157) = {156};
Line Loop(158) = {95, 84, -115};
Plane Surface(159) = {158};
Line Loop(160) = {115, -85, 98};
Plane Surface(161) = {160};
Line Loop(162) = {85, 116, 97};
Plane Surface(163) = {162};
Line Loop(164) = {116, -96, 84};
Plane Surface(165) = {164};
Physical Line(166) = {1, 63, 64, 2, 57, 81};
Physical Line(167) = {7, 8};
Physical Line(168) = {82, 58, 4, 70, 71, 3};
Physical Line(169) = {5, 6};
Physical Surface(1) = {100, 109, 111, 113, 121, 123, 125, 119, 114, 106, 104, 102};
Physical Surface(2) = {127, 133, 131, 129, 137, 139, 141, 135, 143, 149, 147, 145, 153, 155, 157, 151, 159, 165, 163, 161};


EndIf
If(square == 1)
// square
Point(150) = {-sep_norm/2-hx+d/2+radius1, -hy/2+radius1, 0,lc3};
Point(151) = {-sep_norm/2-hx+d/2-radius1, -hy/2+radius1, 0,lc3};
Point(152) = {-sep_norm/2-hx+d/2+radius1, -hy/2-radius1, 0,lc3};
Point(153) = {-sep_norm/2-hx+d/2-radius1, -hy/2-radius1, 0,lc3};
Point(154) = {sep_norm/2-hx+d/2+radius1, -hy/2+radius1, 0,lc3};
Point(155) = {sep_norm/2-hx+d/2-radius1, -hy/2+radius1, 0,lc3};
Point(156) = {sep_norm/2-hx+d/2+radius1, -hy/2-radius1, 0,lc3};
Point(157) = {sep_norm/2-hx+d/2-radius1, -hy/2-radius1, 0,lc3};


Line(21) = {152, 8};
Line(22) = {8, 153};
Line(23) = {153, 7};
Line(24) = {7, 151};
Line(41) = {151, 6};
Line(42) = {6, 150};
Line(43) = {150, 9};
Line(44) = {151, 7};
Line(45) = {7, 153};
Line(46) = {153, 8};
Line(47) = {8, 152};
Line(48) = {152, 9};
Line(49) = {155, 15};
Line(50) = {15, 154};
Line(51) = {154, 18};
Line(52) = {16, 155};
Line(53) = {157, 17};
Line(54) = {17, 156};
Line(55) = {156, 18};
Line(56) = {16, 157};


EndIf


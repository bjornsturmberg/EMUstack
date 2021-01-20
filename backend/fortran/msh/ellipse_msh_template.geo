// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
ff = 0;
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
a1 = 0;
radius1 = (a1/(2*d_in_nm))*d;
ellipticity = 0;
angle = 0;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres
lc4 = lc/1; // centres of top and bottom

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;

// Ellipse aligne with main axes
// so I force the simmetry of the unit cell
If(angle == 0)

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
Line(1) = {1,10};
Line(2) = {10,4};
Line(3) = {2,12};
Line(4) = {12,3};
Line(5) = {1,11};
Line(6) = {11,2};
Line(7) = {4,13};
Line(8) = {13,3};
Line(9) = {11,7};
Line(10) = {7,5};
Line(11) = {5,9};
Line(12) = {9,13};
Line(13) = {10,6};
Line(14) = {6,5};
Line(15) = {5,8};
Line(16) = {8,12};

Ellipsis(17) = {9,5,6,6};
Ellipsis(18) = {6,5,7,7};
Ellipsis(19) = {7,5,8,8};
Ellipsis(20) = {8,5,9,9};

Line Loop(21) = {7,-12,17,-13,2};
Plane Surface(22) = {21};
Line Loop(23) = {1,13,18,-9,-5};
Plane Surface(24) = {23};
Line Loop(25) = {6,3,-16,-19,-9};
Plane Surface(26) = {25};
Line Loop(27) = {4,-8,-12,-20,16};
Plane Surface(28) = {27};
Line Loop(29) = {17,14,11};
Plane Surface(30) = {29};
Line Loop(31) = {18,10,-14};
Plane Surface(32) = {31};
Line Loop(33) = {19,-15,-10};
Plane Surface(34) = {33};
Line Loop(35) = {20,-11,15};
Plane Surface(36) = {35};

Physical Line(1) = {1,2};
Physical Line(2) = {3,4};
Physical Line(3) = {5,6};
Physical Line(4) = {7,8};

Physical Surface(1) = {24,26,28,22};
Physical Surface(2) = {30,32,34,36};

// Symmetry is brocken, so no symmetry forcing
Else

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

Rotate{{0, 0, 1}, {0.5 * d, -0.5 * dy, 0}, angle} {Point{6,7,8,9};}

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Ellipsis(5) = {9,5,6,6};
Ellipsis(6) = {6,5,7,7};
Ellipsis(7) = {7,5,8,8};
Ellipsis(8) = {8,5,9,9};

Line Loop(9) = {1,2,3,4};
Line Loop(10) = {5,6,7,8};
Plane Surface(11) = {9,10};
Plane Surface(12) = {10};

// Add points to surface for mesh control
Point{5} In Surface{12};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Physical Surface(1) = {11};
Physical Surface(2) = {12};

EndIf

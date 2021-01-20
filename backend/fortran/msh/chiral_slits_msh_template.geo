// Template mesh geometry file for a single inclusion.
// By default it will be circular, can also be made to be
// elliptical or square.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;
Geometry.CopyMeshingMethod = 1;

d = 1; // grating period
ff = 0;

// input geo params
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm / d_in_nm;
a1 = 0; // horizontal rectangle side
a2 = 0;  // vertical rectangle side
smooth = 0;
offset_x_nm = 0;
offset_y_nm = 0;
os_x = (offset_x_nm / d_in_nm) * d;
os_y = (offset_y_nm / dy_in_nm) * dy;

// input mesh params
lc = 0;    // 0.501 0.201 0.0701;
lc2 = lc/1;  // on edges
lc3 = lc/1; // on corners

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;

// derived params
r1 = (a1 / (2 * d_in_nm)) * d;
r2 = (a2 / (2 * d_in_nm)) * d;
s2 = r2 * smooth;
r1r = r1 - s2;
r2r = r2 - s2;

// unit cell
Point(1) = {0, 0, 0, lc2}; // first side
Point(2) = {0, -dy, 0, lc2};
Point(3) = {d, -dy, 0, lc2};
Point(4) = {d, 0, 0, lc2};

// rectangle slit
Point(10) = {0.5 * d + r1 + os_x, -0.5 * dy + os_y, 0, lc2}; // upper right
Point(11) = {0.5 * d + r1 + os_x, -0.5 * dy + r2r + os_y, 0, lc3};
Point(12) = {0.5 * d + r1r + os_x, -0.5 * dy + r2r + os_y, 0, lc2};
Point(13) = {0.5 * d + r1r + os_x, -0.5 * dy + r2 + os_y, 0, lc3};
Point(14) = {0.5 * d + os_x, -0.5 * dy + r2 + os_y, 0, lc2};
Point(15) = {0.5 * d - r1r + os_x, -0.5 * dy + r2 + os_y, 0, lc3}; // upper left
Point(16) = {0.5 * d - r1r + os_x, -0.5 * dy + r2r + os_y, 0, lc2};
Point(17) = {0.5 * d - r1 + os_x, -0.5 * dy + r2r + os_y, 0, lc3};
Point(18) = {0.5 * d - r1 + os_x, -0.5 * dy + os_y, 0, lc2};
Point(19) = {0.5 * d - r1 + os_x, -0.5 * dy - r2r + os_y, 0, lc3}; //lower left
Point(20) = {0.5 * d - r1r + os_x, -0.5 * dy - r2r + os_y, 0, lc2};
Point(21) = {0.5 * d - r1r + os_x, -0.5 * dy - r2 + os_y, 0, lc3};
Point(22) = {0.5 * d + os_x, -0.5 * dy - r2 + os_y, 0, lc2};
Point(23) = {0.5 * d + r1r + os_x, -0.5 * dy - r2 + os_y, 0, lc3}; //lower rigth
Point(24) = {0.5 * d + r1r + os_x, -0.5 * dy - r2r + os_y, 0, lc2};
Point(25) = {0.5 * d + r1 + os_x, -0.5 * dy - r2r + os_y, 0, lc3};

// Additional points for mesh control
Point(26) = {0.5 * d, -0.5 * dy, 0, lc};
Point(27) = {0.5 * d + 0.3 * d, -0.5 * dy + 0.3 * dy, 0, lc};
Point(28) = {0.5 * d + 0.3 * d, -0.5 * dy - 0.3 * dy, 0, lc};
Point(29) = {0.5 * d - 0.3 * d, -0.5 * dy + 0.3 * dy, 0, lc};
Point(30) = {0.5 * d - 0.3 * d, -0.5 * dy - 0.3 * dy, 0, lc};

// LINES: unit cell
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// LINES: rectangle slit
Line(21) = {10, 11};
Ellipse(22) = {11, 12, 13, 13};
Line(23) = {13, 14};
Line(24) = {14, 15};
Ellipse(25) = {15, 16, 17, 17};
Line(26) = {17, 18};
Line(27) = {18, 19};
Ellipse(28) = {19, 20, 21, 21};
Line(29) = {21, 22};
Line(30) = {22, 23};
Ellipse(31) = {23, 24, 25, 25};
Line(32) = {25, 10};

// LINES: copies for the other 3 slits
SlitLines2[] = Rotate{{0, 0, 1}, {0.5 * d, -0.5 * dy, 0}, Pi / 2} {Duplicata{Line{21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
}
}
;
SlitLines3[] = Rotate{{0, 0, 1}, {0.5 * d, -0.5 * dy, 0}, Pi} {Duplicata{Line{21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
}
}
;
SlitLines4[] = Rotate{{0, 0, 1}, {0.5 * d, -0.5 * dy, 0}, 3 * Pi / 2} {Duplicata{Line{21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
}
}
;

// LINELOOPS for the unit cell and the four rectangular slits
Curve Loop(41) = {1, 2, 3, 4};
Curve Loop(42) = {21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
Curve Loop(43) = SlitLines2[];
Curve Loop(44) = SlitLines3[];
Curve Loop(45) = SlitLines4[];

// Surfaces for the four slits and the unit cell minus the four slits
Plane Surface(42) = {42};
Plane Surface(43) = {43};
Plane Surface(44) = {44};
Plane Surface(45) = {45};
Plane Surface(41) = {41, 42, 43, 44, 45};

// Add points to surface for mesh control
Point{26} In Surface{41};
Point{27} In Surface{41};
Point{28} In Surface{41};
Point{29} In Surface{41};
Point{30} In Surface{41};

//PHYSICAL ENTITIES
Physical Line(1) = {1}; // external
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Physical Surface(1) = {41}; //internal
Physical Surface(2) = {42, 43, 44, 45};

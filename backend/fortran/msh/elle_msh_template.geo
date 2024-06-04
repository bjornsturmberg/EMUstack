// Template mesh geometry file for ELLE shape.
// See Zanotto et al Nanophotonics 8, 2291 (2019)

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;

// input parameters in units of period d
L1 = 0;
L2 = 0;
L3 = 0;
L4 = 0;


// mesh finesse parameters
lc = 0.100000; // 0.501 0.201 0.0701;
lc2 = lc/2.000000; // mesh at ELLE boundary

// unitary cell
Point(1) = {-0.5*d,-0.5*dy, 0, lc};
Point(2) = {0.5*d,-0.5*dy, 0, lc};
Point(3) = {0.5*d,0.5*dy, 0, lc};
Point(4) = {-0.5*d,0.5*dy, 0, lc};

// L points
Point(5)  = {-0.5*L1,    -0.5*L2,    0, lc2}; // lower left corner
Point(6)  = { 0.5*L1,    -0.5*L2,    0, lc2}; // lower right corner
Point(7)  = { 0.5*L1,    -0.5*L2+L4, 0, lc2};
Point(8)  = {-0.5*L1+L3, -0.5*L2+L4, 0, lc2};
Point(9)  = {-0.5*L1+L3,  0.5*L2,    0, lc2};
Point(10) = {-0.5*L1,     0.5*L2,    0, lc2}; // upper left corner


// L lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};
Line(8)  = {8,9};
Line(9)  = {9,10};
Line(10) = {10,5};


// Line loops defining the surfaces
Line Loop(1) = {5,6,7,8,9,10};
Line Loop(2) = {1,2,3,4};

// Surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2,1};

// Physical Surfaces where materials get assigned
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Physical Surface(1) = {2}; //internal
Physical Surface(2) = {1};
//+
Show "*";
//+
Show "*";

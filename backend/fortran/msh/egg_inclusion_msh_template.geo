// Template mesh geometry file egg plus active material inclusion
// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
ff = 0;
d_in_nm = 0;
dy_in_nm = 0;
dy = dy_in_nm/d_in_nm;
r1 = 0;
r2 = 0;
r3 = 0;
angle = 0;
shiftx = 0;
shifty = 0;
///////////////////////////////////
// active material inclusion params
///////////////////////////////////
r_inc = 0;
inc_shiftx = 0;
inc_shifty = 0;
///////////////////////////////////
///////////////////////////////////
///////////////////////////////////
rr1 = (r1/(d_in_nm))*d;
rr2 = (r2/(d_in_nm))*d;
rr3 = (r3/(d_in_nm))*d;
rr_inc = (r_inc/(d_in_nm))*d;
inc_sshiftx = (inc_shiftx/(d_in_nm))*d;
inc_sshifty = (inc_shifty/(d_in_nm))*d;
inc_delta = 0.0000001; // this is a numerical trick to distinguish point(5) from point(105) when there is no shift
sshiftx=(shiftx/(d_in_nm))*d;
sshifty=(shifty/(d_in_nm))*d;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres
///////////////////////////////////
// active material inclusion params
///////////////////////////////////
lc4 = lc/1; // active material inclusion
///////////////////////////////////
///////////////////////////////////
///////////////////////////////////

hy = dy; // Thickness: Square profile => hy=d
hx = 0.;

// Simplified unit cell
Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// circular inclusion of active material, I include sshiftx and sshifty so it moves together with the egg
Point(105) = {inc_sshiftx + sshiftx - hx+d/2 + inc_delta,inc_sshifty +  -hy/2 + sshifty, 0,lc4};
Point(106) = {inc_sshiftx + sshiftx - hx+d/2,inc_sshifty + -hy/2 - rr_inc + sshifty, 0, lc4};
Point(107) = {inc_sshiftx + sshiftx - hx+d/2 - rr_inc,inc_sshifty + -hy/2 + sshifty, 0, lc4};
Point(108) = {inc_sshiftx + sshiftx - hx+d/2,inc_sshifty + -hy/2 + rr_inc + sshifty, 0, lc4};
Point(109) = {inc_sshiftx + sshiftx - hx+d/2 + rr_inc,inc_sshifty + -hy/2 + sshifty, 0, lc4};

// Vertices of the egg
Point(5) = {sshiftx - hx+d/2, -hy/2 + sshifty, 0,lc3};
Point(6) = {sshiftx - hx+d/2, -hy/2 - rr3 + sshifty, 0, lc2};
Point(7) = {sshiftx - hx+d/2 - rr1, -hy/2 + sshifty, 0, lc2};
Point(8) = {sshiftx - hx+d/2, -hy/2 + rr2 + sshifty, 0, lc2};
Point(9) = {sshiftx - hx+d/2 + rr1, -hy/2 + sshifty, 0, lc2};

Rotate{{0, 0, 1}, {sshiftx - hx+d/2, -hy/2 + sshifty, 0}, angle} {Point{6,7,8,9,105,106,107,108,109};}

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Ellipsis(5) = {9,5,6,6};
Ellipsis(6) = {6,5,7,7};
Ellipsis(7) = {7,5,8,8};
Ellipsis(8) = {8,5,9,9};

// active materials inclusion arcs
Ellipsis(105) = {109,105,106,106};
Ellipsis(106) = {106,105,107,107};
Ellipsis(107) = {107,105,108,108};
Ellipsis(108) = {108,105,109,109};


Line Loop(9) = {1,2,3,4};
Line Loop(10) = {5,6,7,8};
Line Loop(100) = {105,106,107,108};
Plane Surface(11) = {9,10};
Plane Surface(12) = {10,100};
Plane Surface(101) = {100};

// Add points to surface for mesh control
Point{5} In Surface{12};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Physical Surface(1) = {11};
Physical Surface(2) = {12};
Physical Surface(3) = {101};
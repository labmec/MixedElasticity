//+
//SetFactory("OpenCASCADE");
a=0.2;
Point(1) = {-1,-0.5,0,a};
Point(2) = { 1,-0.5,0,a};
Point(3) = { 1, 0.5,0,a};
Point(4) = {-1, 0.5,0,a};

Point(10) = {-0.05,0,0,a*0.05};
Point(11) = {0.05,0,0,a*0.05};
Point(12) = {0.0,0,0,a};
//+
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {10,12, 11};
Circle(6) = {11,12, 10};

//+
Curve Loop(1) = {4, 1, 2, 3,-5,-6};
//+
Plane Surface(1) = {1};
Recombine Surface{1};
//+
Physical Surface("Domain") = {1};
//+
Physical Curve("Lower") = {1};

//+
Physical Curve("Right") = {2};
Physical Curve("Upper") = {3};
Physical Curve("Left") = {4};
//+
Physical Curve("Hole") = {5,6};
//+
Physical Point("Corner") = {1};

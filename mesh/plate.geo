//+
//SetFactory("OpenCASCADE");
a=0.2;
Point(1) = {-1,-0.25,0,a};
Point(2) = { 1,-0.25,0,a};
Point(3) = { 1, 0.25,0,a};
Point(4) = {-1, 0.25,0,a};

//+
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
n = 0;
Transfinite Curve {1, 3} = 4*n+1 Using Progression 1;
Transfinite Curve {2, 4} = n+1 Using Progression 1;
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
Transfinite Surface{1};
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
Physical Point("Corner") = {1};
//+

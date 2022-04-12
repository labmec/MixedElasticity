//+
//SetFactory("OpenCASCADE");
a=0.1;
Point(1) = {0,0,0,a};
Point(2) = {1,0,0,a};
Point(3) = {1,2,0,a};
Point(4) = {0,2,0,a};

Point(10) = {0.45,1,0,a};
Point(11) = {0.55,1,0,a};
//+
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {10, 11};
Line(6) = {11, 10};

//+
Curve Loop(1) = {4, 1, 2, 3,5,6};
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
Physical Curve("Fracture") = {5,6};
//+
Physical Point("FractureCorner") = {10, 11};

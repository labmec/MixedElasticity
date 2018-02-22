//+
SetFactory("OpenCASCADE");
angle = 2* Pi/9;
cangle = Pi/2-angle;
angle2 = Pi/6;
cangle2 = Pi/2-angle2;
Rc = 15;
thick = 0.06;
Rd = Rc/Sin(angle);
R1 = Rd-thick/2;
R2 = Rd+thick/2;
Point(0) = {0, 0, 0, 1.0};
Point(1) = {0, R1, 0, 0.5};
Point(2) = {0, R2, 0, 0.5};
msh=0.1;
msh2 = 0.5;
mshverysmall=0.01;
Point(3) = {R1*Cos(cangle), R1 * Sin(cangle), 0, mshverysmall};
Point(4) = {R2 * Cos(cangle), R2 * Sin(cangle), 0, mshverysmall};
Point(5) = {R1*Cos(cangle) + 0.6, R2 * Sin(cangle) , 0, msh2};
Point(6) = {R1*Cos(cangle), R2 * Sin(cangle) - 0.5, 0, msh2};
Point(7) = {R1*Cos(cangle) + 0.6, R2 * Sin(cangle) - 0.5,0, msh2};
Point(8) = {R1 * Cos(cangle2), R1 * Sin(cangle2), 0, msh};
Point(9) = {R2 * Cos(cangle2), R2 * Sin(cangle2), 0, msh};

//+
Circle(1) = {8, 0, 1};
//+
Circle(2) = {9, 0, 2};
//+
Circle(3) = {3, 0, 8};
//+
Circle(4) = {4, 0, 9};
//+
//+
Line(5) = {2, 1};
//+
Line(6) = {3, 4};
//+
Line(7) = {3, 6};
//+
Line(8) = {6, 7};
//+
Line(9) = {7, 5};
//+
Line(10) = {5, 4};
//+
Line(11) = {9, 8};
//+
Line Loop(1) = {7, 8, 9, 10, -6};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {3, 1, -5, -2, -4, -6};
//+
Plane Surface(2) = {2};
//+
Physical Surface("ELASTICITY1") = {2};
//+
Physical Surface("ELASTICITY2") = {1};
//+
Physical Line("SUPPORT") = {8};
//+
Physical Line("ZERO") = {9, 10, 7, 3, 4, 2, 1};
//+
Physical Line("MOMENT") = {6};

Physical Line("SYMMETRY") = {5};

Physical Point("POINTBC") = {6};
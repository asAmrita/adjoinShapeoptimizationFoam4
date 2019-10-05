Lx = 1;
Ly = 1;

//+
Point(1) = {0, 0, 0, 0};
//+
Point(2) = {Lx, 0, 0, 0};
//+
Point(3) = {Lx, Ly, 0, 0};
//+
Point(4) = {0, Ly, 0, 0};
//+
Point(5) = {0.5, 0, 0, 0};
//+
Point(6) = {0.75, 0, 0, 0};
//+
Point(7) = {Lx, 0.5, 0, 0};
//+
Point(8) = {Lx, 0.75, 0, 0};
//+
Line(1) = {1, 5};
//+
Line(2) = {5, 6};
//+
Line(3) = {6, 2};
//+
Line(4) = {2, 7};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 1};
//+
Line Loop(1) = {8, 1, 2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};

//+
Transfinite Line {8} = 81 Using Progression 1;
//+
Transfinite Line {4} = 41 Using Progression 1;
//+
Transfinite Line {5} = 21 Using Progression 1;
//+
Transfinite Line {6} = 21 Using Progression 1;
//+
Transfinite Line {3} = 21 Using Progression 1;
//+
Transfinite Line {2} = 21 Using Progression 1;
//+
Transfinite Line {1} = 41 Using Progression 1;
//+
Transfinite Line {7} = 81 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {1} = {4, 1, 2, 3};
//+
Recombine Surface {1};
//+
Extrude {0, 0, 0.1} {
  Surface{1}; Layers{1}; Recombine;
}
//+
Physical Surface("inlet") = {29};
//+
Physical Surface("outlet") = {41};
//+
Physical Surface("upperWall") = {21,49,45,37};

//+
Physical Surface("lowerWall") = {45,37,25,33};

Physical Surface("frontAndBack") = {50, 1};
//+
Physical Volume("body") = {1};

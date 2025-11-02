Mesh.CharacteristicLengthMin = 3.8;
Mesh.CharacteristicLengthMax = 3.9;

Merge "corpo-peri.stl";
Merge "corpo-endo.stl";
Merge "colo-peri.stl";
Merge "colo-endo.stl";
Merge "colo-base.stl";
Merge "corpo-base.stl";
//Merge "colo-topo.stl";


Coherence Mesh;

CreateTopology;

CreateGeometry;


//+
Physical Surface("peri-corpo", 30) = {1};
//+
Physical Surface("endo-corpo", 20) = {2};
//+
Physical Surface("base-corpo", 10) = {6};
//+
Physical Surface("peri-colo", 60) = {3};
//+
Physical Surface("endo-colo", 40) = {4};
//+
Physical Surface("base-colo", 50) = {5};
//+
Surface Loop(1) = {1, 2, 6};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {3, 4, 5, -6};
//+
Volume(2) = {2};
//+
Physical Volume("corpo", 1) = {1};
//+
Physical Volume("colo", 2) = {2};




// Conventional cell (no 3D current collector) matching NEW mesh dimensions
// Domain: 256.401 x 256.401 x 228 µm
// Cathode: z = 0 to 197.757 µm
// Separator: z = 197.757 to 228 µm
// cat_cc: flat plane at z = 0 (bottom face)
// top: flat plane at z = 228 (top face)

// Mesh size — adjust for desired refinement
lc = 20.0;  // characteristic element size in µm

// Dimensions
Lx = 256.401;
Ly = 256.401;
Lz_cathode = 197.757;   // cathode thickness (matches new mesh separator start)
Lz_total = 228;         // total cell thickness

// --- Bottom face (z=0, will be cat_cc) ---
Point(1) = {0, 0, 0, lc};
Point(2) = {Lx, 0, 0, lc};
Point(3) = {Lx, Ly, 0, lc};
Point(4) = {0, Ly, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- Cathode volume: extrude bottom face to z = Lz_cathode ---
cathode_extrude[] = Extrude {0, 0, Lz_cathode} { Surface{1}; };
// cathode_extrude[0] = top surface (cathode/separator interface)
// cathode_extrude[1] = volume

// --- Separator volume: extrude interface to z = Lz_total ---
separator_extrude[] = Extrude {0, 0, Lz_total - Lz_cathode} { Surface{cathode_extrude[0]}; };
// separator_extrude[0] = top surface (top boundary)
// separator_extrude[1] = volume

// --- Physical groups ---
// Volumes
Physical Volume("cathode") = {cathode_extrude[1]};
Physical Volume("block_0") = {separator_extrude[1]};

// Boundaries
Physical Surface("cat_cc") = {1};                        // bottom face, z=0
Physical Surface("top") = {separator_extrude[0]};        // top face, z=228

// Generate 3D mesh with tetrahedra
Mesh.Algorithm3D = 1;  // Delaunay
Mesh 3;

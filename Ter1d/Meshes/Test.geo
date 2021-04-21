// Gmsh project created on Tue Apr 20 21:39:48 2021
h = 0.5;
L=5;
Point (5) = {0,L,0,h};
Point (6) = {0,-L,0,h};
Point (7) = {0,0,0,h};
Point (1) = {-L,-L,0,h};
Point (2) = {-L,L,0,h};
Point (3) = {L,L,0,h};
Point (4) = {L,-L,0,h};

Line (1) = {1 ,2}; // bas
Line (2) = {2 ,3}; // droite
Line (3) = {3 ,4}; // haut
Line (4) = {4 ,1}; // gauche

Line Loop (1) = {1 ,2 ,3 ,4};

Plane Surface(1)={1};
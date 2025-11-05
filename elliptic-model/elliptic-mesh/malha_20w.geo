// Parametros do modelo
xc = 0.0;
yc = 0.0;
zc = 0.0;
a = 123.0/2.0;
b = 135.0/2.0;
c = 135.0/2.0;
thick = 8.4;
m_cut = 15.0;

// Tamanho malha
msize = 8.0;

// Center
Point(1) = {xc, yc, zc, msize};
// Pontos elipse interna plano xy
Point(2) = {a, yc, zc, msize};
Point(3) = {-a, yc, zc, msize};
Point(4) = {xc, b, zc, msize};
Point(5) = {xc, -b, zc, msize};
// Pontos elipse externa plano xy
Point(6) = {a+thick, yc, zc, msize};
Point(7) = {-a-thick, yc, zc, msize};
Point(8) = {xc, b+thick, zc, msize};
Point(9) = {xc, -b-thick, zc, msize};

// Arcos de Elipse:
// Ellipse(n) = {p_inicial, p_central, p_eixoPrincipal,p_final}
/////////////////////////////////////////////////////////////
// Fundo do utero
// Elipse interna plano xy
Ellipse(1) = {2, 1, 4, 4};
Ellipse(2) = {4, 1, 3, 3};
Ellipse(3) = {3, 1, 5, 5};
Ellipse(4) = {5, 1, 2, 2};

// Elipse externa plano xy
Ellipse(5) = {6, 1, 8, 8};
Ellipse(6) = {8, 1, 7, 7};
Ellipse(7) = {7, 1, 9, 9};
Ellipse(8) = {9, 1, 6, 6};

// Ponto elipse interna plano xz
Point(10) = {xc, yc, c, msize};
// Ponto elipse externa plano xz
Point(11) = {xc, yc, c+thick, msize};

// Elipses interna plano xz
Ellipse(9) = {2, 1, 10, 10};
Ellipse(10) = {10, 1, 3, 3};

// Elipse externa plano xz
Ellipse(11) = {6, 1, 11, 11};
Ellipse(12) = {11, 1, 7, 7};

// Elipses interna plano yz
Ellipse(13) = {4, 1, 10, 10};
Ellipse(14) = {10, 1, 5, 5};
// Elipse externa plano yz
Ellipse(15) = {8, 1, 11, 11};
Ellipse(16) = {11, 1, 9, 9};

//////////////////////////////////////////////////////////
// Parte inferior do utero
// Inferior interna plano xz
xb_i = 38.6551707048645980;
Point(12) = {xb_i,yc,-c+m_cut,msize};
Point(13) = {-xb_i,yc,-c+m_cut,msize};
Ellipse(17) = {12,1,2,2};
Ellipse(18) = {3,1,3,13};
// Inferior interna plano yz
yb_i = 42.42640687119285146;
Point(14) = {xc,yb_i,-c+m_cut,msize};
Point(15) = {xc,-yb_i,-c+m_cut,msize};
Ellipse(19) = {14,1,4,4};
Ellipse(20) = {5,1,15,15};

// Inferior externa plano xz
xb_e = 50.48075485403204379;
Point(16) = {xb_e,yc,-c+m_cut,msize};
Point(17) = {-xb_e,yc,-c+m_cut,msize};
Ellipse(21) = {16,1,6,6};
Ellipse(22) = {7,1,7,17};
// Inferior externa plano yz
yb_e = 54.81386685866998746;
Point(18) = {xc,yb_e,-c+m_cut,msize};
Point(19) = {xc,-yb_e,-c+m_cut,msize};
Ellipse(23) = {18,1,8,8};
Ellipse(24) = {9,1,19,19};

//////////////////////////////////////
// Base
// Centro base
Point(20) = {xc,yc,-c+m_cut,msize};
//Base interna
Ellipse(25) = {12,20,14,14};
Ellipse(26) = {14,20,13,13};
Ellipse(27) = {13,20,15,15};
Ellipse(28) = {15,20,12,12};
// Base externa
Ellipse(29) = {16,20,18,18};
Ellipse(30) = {18,20,17,17};
Ellipse(31) = {17,20,19,19};
Ellipse(32) = {19,20,16,16};

// Definindo Superficies
// curve loop(n) = {rotulo das fronteiras da superficie}
// Endometrio 1o octante:
Curve Loop(1) = {1,13,-9}; //interno
Surface(1) = {1};
Curve Loop(2) = {5,15,-11}; //externo
Surface(2) = {2};

// Endometrio 2o octante:
Curve Loop(3) = {-2,13,10}; //interno
Surface(3) = {3};
Curve Loop(4) = {-6,15,12}; //externo
Surface(4) = {4};

// Endometrio 3o octante:
Curve Loop(5) = {3,-14,10}; //interno
Surface(5) = {5};
Curve Loop(6) = {7,-16,12}; //externo
Surface(6) = {6};

// Endometrio 4o octante:
Curve Loop(7) = {-4,-14,-9}; //interno
Surface(7) = {7};
Curve Loop(8) = {-8,-16,-11}; //externo
Surface(8) = {8};

// Endometrio 5o octante:
Curve Loop(9) = {25,19,-1,-17}; //interno
Surface(9) = {9};
Curve Loop(10) = {29,23,-5,-21}; //externo
Surface(10) = {10};

// Endometrio 6o octante:
Curve Loop(11) = {-26,19,2,18}; //interno
Surface(11) = {11};
Curve Loop(12) = {-30,23,6,22}; //externo
Surface(12) = {12};

// Endometrio 7o octante:
Curve Loop(13) = {27,-20,-3,18}; //interno
Surface(13) = {13};
Curve Loop(14) = {31,-24,-7,22}; //externo
Surface(14) = {14};

// Endometrio 8o octante:
Curve Loop(15) = {-28,-20,4,-17}; //interno
Surface(15) = {15};
Curve Loop(16) = {-32,-24,8,-21}; //externo
Surface(16) = {16};

// // Base em z=-c+cut
Curve Loop(17) = {-28,-27,-26,-25};
Curve Loop(18) = {-32,-31,-30,-29};
Plane Surface(17) = {17,18};

// Cria superficie do corpo
Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};

Volume(1) = {1};

// Marca as surperficies com os rotulos necessarios
Physical Surface("Base", 10) = {17};

Physical Surface("Endometrio", 20) = {1, 3, 5, 7, 9, 11, 13, 15};

Physical Surface("Perimetrio", 30) = {2, 4, 6, 8, 10, 12, 14, 16};

Physical Volume("Vol") = {1};
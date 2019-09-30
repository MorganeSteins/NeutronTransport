#include "aleat_tools.hpp"
#include <math.h>
#include<iostream>
#include <vector>
#include <cstdlib>
#include <typeinfo>

using namespace std;

//tirage d'une variable aléatoire pour le parcourt entre x_n-1 et x_n
double parcourt_x(point p) {
    double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); // loi uniforme sur (0,1)
    return -(p.get_mu()/p.sigmaT)*log(1-y);
}

//nouveau point p après déplacement
point deplacement_x(point p){
    double new_x =  p.get_x()+p.get_mu()*parcourt_x(p);
    point p2;
    p2.set_mu(p.get_mu());
    p2.set_x(new_x);
    return p2;
}

//tirage d'une variable aléatoire pour la continuation de la marche 0 arret, 1 continuer
double do_I_stop(double sigmaS, double sigmaA) {
    double y = (static_cast <double> (rand()) / static_cast <double> (RAND_MAX))*(sigmaA+sigmaS); 
    if (y<sigmaA) {return 0;}
    else {return 1;}
}

//tirage d'un cos(angle) uniforme entre -1 et 1
double new_mu() {
    return 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX))-1.;
}

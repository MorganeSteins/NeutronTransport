#include "MC.hpp"
#include <math.h>
#include<iostream>

using namespace std;



double parcourt_x(point p, double sigmat) {
    return -(p.get_mu()/sigmat)*log(1-p.get_x());
}

double new_x(point p, double sigmat){
    double deplacmt = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    return x+parcourt_x(p, sigma);
}


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

list_points::list_points(list_points &list_to_copy) {
    dim = list_to_copy.dim;
    val = new point(dim);
    for (int i=0;i<V.dim;i++) {this->points[i]=list_to_copy.points[i];}

}

list_points::~list_points() {
    delete points;
    dim = 0;
}

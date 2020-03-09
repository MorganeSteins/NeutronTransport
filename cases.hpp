#ifndef CASES_HPP_INCLUDED
#define CASES_HPP_INCLUDED

#include <vector>
#include <math.h>
#include <iostream>
#include "points.hpp"
#include "aleat_tools.hpp"

class source_point {
    public :
    point change_x(point p) ;
} ;

point source_point::change_x(point p){
    return parcours_x_p(point(0,p.get_mu()));
}

class source_uniforme {
    public :
    point change_x(point p) ;
};

point source_uniforme::change_x(point p){
    return deplacement_x(point(new_x(),p.get_mu()));
}

template <typename SOURCE > 
double density_no_scattering_MC(double mu, double x, int N)
{
    double dx = 1. / 100;
    double rsl = 0.;
    double pt = 0.;
    SOURCE source ;
    for (int i = 0; i < N; i++)
    {
        pt = source.change_x(point(x,mu)).get_x();
        if (pt > x - dx / 2 && pt < dx / 2 + x) rsl++;
    }
    return rsl / (N * point::sigmaT * dx);
}

#endif // CASES_HPP_INCLUDED
#ifndef ALEAT_TOOLS_HPP_INCLUDED
#define ALEAT_TOOLS_HPP_INCLUDED

#include <vector>
#include <iostream>
#include "points.hpp"

using namespace std;


double parcourt_x(point p, double sigmat);
point deplacement_x(point p, double sigmat);
double do_I_stop(double sigmaS, double sigmaA);
double new_mu() ;
vector_points no_absorption_MC(int N, double mu);


#endif // ALEAT_TOOLS_HPP_INCLUDED
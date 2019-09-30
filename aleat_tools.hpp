#ifndef ALEAT_TOOLS_HPP_INCLUDED
#define ALEAT_TOOLS_HPP_INCLUDED

#include <vector>
#include <iostream>

using namespace std;

void print_vector(vector_points &v);
double parcourt_x(point p, double sigmat);
point deplacement_x(point p, double sigmat);
double do_I_stop(double sigmaS, double sigmaA);
double new_mu() ;


#endif // ALEAT_TOOLS_HPP_INCLUDED
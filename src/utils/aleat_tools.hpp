#ifndef ALEAT_TOOLS_HPP_INCLUDED
#define ALEAT_TOOLS_HPP_INCLUDED

#include <vector>
#include <iostream>
#include "points.hpp"

using namespace std;

double new_x();
double parcours_x(point p);
point parcours_x_p(point p);
point deplacement_x(point p);
int do_I_stop(double sigmaS, double sigmaA);
int do_I_stop_woodcock(double sigmaTmax, double sigmaS, double sigmaT);
double new_mu() ;

double vector_max(vector<double> v);
double compare_vect(vector<double> v1, vector<double> v2);

#endif // ALEAT_TOOLS_HPP_INCLUDED
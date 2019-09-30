#ifndef MC_HPP_INCLUDED
#define MC_HPP_INCLUDED

#include <vector>
#include <math.h>
#include <iostream>
#include "points.hpp"
#include "aleat_tools.hpp"

using namespace std;

vector_points no_scattering_MC(int N, double mu);

double density_no_scattering_MC(point p, int N);

double density_no_scattering(point p); 

#endif // MC_HPP_INCLUDED
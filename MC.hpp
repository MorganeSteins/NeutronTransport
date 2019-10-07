#ifndef MC_HPP_INCLUDED
#define MC_HPP_INCLUDED

#include <vector>
#include <math.h>
#include <iostream>
#include "points.hpp"
#include "aleat_tools.hpp"

using namespace std;

/* Fonctions pour le cas homogène et source ponctuelle*/

vector_points no_scattering_homog_point_MC(int N, double mu);

double density_no_scattering_homog_point_MC(point p, int N);

double density_no_scattering_homog_point(point p); 



/* Fonctions pour le cas homogène et source uniforme */

vector_points no_scattering_homog_point_MC(int N, point p);

double density_no_scattering_homog_unif_MC(point p, int N);

double density_no_scattering_homog_unif(point p);

#endif // MC_HPP_INCLUDED
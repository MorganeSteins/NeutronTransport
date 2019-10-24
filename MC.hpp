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

vector<double> density_segment_no_scattering_homog_point_MC(int N, int nb_points, double mu);

double density_tilda_no_scattering_homog_point_MC(point p, int N);

double density_no_scattering_homog_point(point p); 



/* Fonctions pour le cas homogène et source uniforme */

vector_points no_scattering_homog_point_MC(int N, point p);

double density_no_scattering_homog_unif_MC(point p, int N);

double density_no_scattering_homog_unif(point p);

vector<double> density_segment_no_scattering_homog_unif_MC(int N, int nb_points, point p);



/* Fonctions pour le cas homogène avec scattering et source uniforme*/

double density_scattering_homog_unif_MC(int N, point p, int max_iter=100000, double epsilon = 0.0001);

vector<double> density_tilda_segment_scattering_homg_unif_MC(int N, int nb_points,int max_iter=100000, double epsilon = 0.0001);

/* Fonctions pour le cas Woodcock et source uniforme */

vector<double> density_tilda_segment_scattering_woodcock_unif_MC(int N, int nb_intervalles, int max_iter, double epsilon );

#endif // MC_HPP_INCLUDED
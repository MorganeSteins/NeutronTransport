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



/* Fonctions pour le cas homogène avec scattering et source uniforme*/

void move_scattering(vector_points start, vector_points final);

/* Densité moyennée en mu au point x*/
double density_tilda_scattering_homg_unif_MC(int N, double x, int max_iter=100000, double epsilon = 0.01);



void save_error(int N_max,point p);
void save_points(int N,int nb_points,double mu, string filename_);


#endif // MC_HPP_INCLUDED
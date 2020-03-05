#ifndef MC_ITER_HPP_INCLUDED
#define MC_ITER_HPP_INCLUDED

#include "points.hpp"
#include "aleat_tools.hpp"

using namespace std ;

/* Objectives : Faster MC implementation with iteration instead of vectors*/
double density_no_scattering_homog_point_MC_iter(point p, int N) ; 
double density_tilda_no_scattering_homog_point_MC_iter(point p, int N) ;
double density_tilda_no_scattering_homog_point_MC_iter(double x, int N) ;
double density_no_scattering_homog_unif_MC_iter(double mu, double x, int N) ;

vector<double> density_segment_no_scattering_homog_unif_MC_iter(int N, int nb_points, double x, double mu) ;


#endif // MC_ITER_HPP_INCLUDED
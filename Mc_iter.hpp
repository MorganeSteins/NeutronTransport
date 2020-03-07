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

point step_homog_unif(double x, double mu, int max_iter=1000) ;
vector<double> density_segment_no_scattering_homog_unif_MC_iter(int N, int nb_points, double x, double mu) ;

vector<double> density_tilda_segment_scattering_homg_unif_MC_iter(int N, int nb_intervalles, int max_iter) ;
vector<double> density_tilda_segment_scattering_woodcock_unif_MC_iter(int N, int nb_intervalles, int max_iter) ;
vector<double> density_tilda_segment_scattering_woodcock_point_MC_iter(int N, int nb_intervalles, int max_iter) ;





#endif // MC_ITER_HPP_INCLUDED
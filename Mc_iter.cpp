#include "Mc_iter.hpp"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std ;

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source ponctuelle isotrope en 0 et matériau homogène.
double density_no_scattering_homog_point_MC_iter(point p, int N)
{
    double dx = 1. / 100;
    double rsl = 0.;
    double pt = 0.;
    for (int i = 0; i < N; i++) // tirer un point et regarder s'il est dans l'intervalle
    {
        pt = parcours_x(point(0, p.get_mu()));
        if (pt > p.get_x() - dx / 2 && pt < dx / 2 + p.get_x()) rsl++;
    }
    return rsl / (N * p.sigmaT * dx);
}

// calcul de la densité neutronique au point x moyennée en mu avec N tirages
// avec une source ponctuelle isotrope en 0 et matériau homogène.
double density_tilda_no_scattering_homog_point_MC_iter(point p, int N)
{
    double dx = 1. / 100;
    double rsl = 0.;
    double pt = 0.;
    for (int i = 0; i < N; i++) // tirer un point et regarder s'il est dans l'intervalle
    {
        pt = parcours_x(point(0, new_mu()));
        if (pt > p.get_x() - dx / 2 && pt < dx / 2 + p.get_x()) rsl++; 
    }
    return rsl / (N * p.sigmaT * dx);
}

double density_tilda_no_scattering_homog_point_MC_iter(double x, int N)
{
    double dx = 1. / 100;
    double rsl = 0.;
    double pt = 0.;
    for (int i = 0; i < N; i++) // tirer un point et regarder s'il est dans l'intervalle
    {
        pt = parcours_x(point(0, new_mu()));
        if (pt > x - dx / 2 && pt < dx / 2 + x) rsl++;
    }
    return rsl / (N * point::sigmaT * dx);
}

double density_no_scattering_homog_unif_MC_iter(double mu, double x, int N)
{
    double dx = 1. / 100;
    double rsl = 0.;
    double pt = 0.;
    for (int i = 0; i < N; i++)
    {
        pt = deplacement_x(point(new_x(), mu)).get_x();
        if (pt > x - dx / 2 && pt < dx / 2 + x) rsl++;
    }
    return rsl / (N * point::sigmaT * dx);
}

vector<double> density_segment_no_scattering_homog_unif_MC_iter(int N, int nb_points, double x, double mu)
{
    double dx = 1. / nb_points;         //largeur des intervalles
    vector<double> freq(nb_points, 0.); //stocke la fréqence pour chaque intervalle
    int indice = 0;
    double pt = 0.; //initialisation

    for (int i = 0; i < N; i++)
    {
        pt = deplacement_x(point(new_x(), mu)).get_x();
        indice = floor((pt / dx)); //indice de l'intervalle où est le neutron
        if (indice < nb_points)
        {
            freq[indice] += (1. / (point::sigmaT * dx)) * (1. / (float)N);
        }
    }
    return freq;
}

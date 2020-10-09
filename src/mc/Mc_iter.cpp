#include "Mc_iter.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>

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
        if (indice < nb_points) freq[indice] += (1. / (point::sigmaT * dx)) * (1. / (float)N);
    }
    return freq;
}

vector<double> density_tilda_segment_scattering_homg_unif_MC_iter(int N, int nb_intervalles, int max_iter)
{
    double dx = 1. / nb_intervalles;
    vector<double> freq(nb_intervalles, 0.);

    //initialisations des variables
    int step = 0;
    int indice = 0;
    point new_point ;
    int stop = 1 ;
    int max_ = 0 ;

    for (int i=0;i<N;i++){
        new_point = point(new_x(), new_mu()) ;
        
        while (step < max_iter && stop != 0){
            new_point = point(deplacement_x(new_point).get_x(), new_mu());
            if (new_point.get_x() >= 0 && new_point.get_x() <= 1 && do_I_stop(point::sigmaS, point::sigmaA)==1) 
            {
                indice = floor(new_point.get_x() / dx); 
                freq[indice] += (1. / (point::sigmaT * dx)) * (1. / (float)N);
            }
            else stop = 0 ;
            step ++ ;
        }
        max_ = max(max_,step) ;
        step = 0 ;  
        stop = 1 ; 
        
    }
    cout<<"Nombre max d'itérations "<<max_<<"\n" ;
    return freq;
}

vector<double> density_tilda_segment_scattering_woodcock_unif_MC_iter(int N, int nb_intervalles, int max_iter)
{
    double dx = 1. / nb_intervalles;
    vector<double> freq(nb_intervalles, 0.);

    //initialisations des variables
    int step = 0;
    int indice = 0;
    point new_point ;
    int stop = 1 ;
    int max_ = 0 ;
    double SigmaT_loc ;

    for (int i=0;i<N;i++){
        new_point = point(new_x(), new_mu()) ;
        
        while (step < max_iter && stop != 0){
            SigmaT_loc = sigmaT_non_cst(new_point) ;
            stop = do_I_stop_woodcock(point::sigmaT, point::sigmaS, SigmaT_loc) ;
            if (stop==2) new_point = deplacement_x(new_point);
            if (stop==1) new_point = point(deplacement_x(new_point).get_x(), new_mu());
            if (new_point.get_x() >= 0 && new_point.get_x() <= 1 && stop==1) 
            {
                indice = floor(new_point.get_x() / dx); 
                freq[indice] += (1. / (SigmaT_loc* dx)) * (1. / (float)N);
            }
            else stop = 0 ;
            step ++ ;
        }
        max_ = max(max_,step) ;
        step = 0 ;  
        stop = 1 ; 
        
    }
    cout<<"Nombre max d'itérations "<<max_<<"\n" ;
    return freq;
}

vector<double> density_tilda_segment_scattering_woodcock_point_MC_iter(int N, int nb_intervalles, int max_iter)
{
    double dx = 1. / nb_intervalles;
    // vector<double> phi(nb_intervalles,0.); // valeurs de l'estimateur
    vector<double> freq(nb_intervalles, 0.);

    //initialisations des variables
    int step = 0;
    int indice = 0;
    point new_point ;
    int stop = 1 ;
    int max_ = 0 ;
    double SigmaT_loc ;

    for (int i=0;i<N;i++){
        new_point = point(0, new_mu()) ;
        
        while (step < max_iter && stop != 0){
            SigmaT_loc = sigmaT_non_cst(new_point) ;
            stop = do_I_stop_woodcock(point::sigmaT, point::sigmaS, SigmaT_loc) ;
            if (stop==2) new_point = deplacement_x(new_point);
            if (stop==1) new_point = point(deplacement_x(new_point).get_x(), new_mu());
            if (new_point.get_x() >= 0 && new_point.get_x() <= 1 && stop==1) 
            {
                indice = floor(new_point.get_x() / dx); 
                freq[indice] += (1. / (SigmaT_loc* dx)) * (1. / (float)N);
            }
            else stop = 0 ;
            step ++ ;
        }
        max_ = max(max_,step) ;
        step = 0 ;  
        stop = 1 ; 
        
    }
    cout<<"Nombre max d'itérations "<<max_<<"\n" ;
    return freq;
}
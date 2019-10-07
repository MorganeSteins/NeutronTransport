#include "MC.hpp"
#include <math.h>
#include<iostream>
#include <vector>
#include <cstdlib>

using namespace std;

/* CAS HOMOGENE SOURCE PONCTUELLE */
//tirage de départ de neutrons dans un matériau homogène, source ponctuelle en 0.
vector_points no_scattering_homog_point_MC(int N, double mu){
    vector_points selection(N);
    for (int i=0; i<N; i++){
        selection.points[i] = point(parcours_x(point(0,mu)),mu);
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source ponctuelle isotrope en 0 et matériau homogène.
double density_no_scattering_homog_point_MC(point p, int N){
    double dx=1./100;
    vector_points selection(N) ;
    selection = no_scattering_homog_point_MC(N,p.get_mu());
    double rsl = 0.;
    for (int i=0;i<N;i++){
        if (selection.points[i].get_x()>p.get_x() && selection.points[i].get_x()<dx+p.get_x()) {
            rsl++;
        }
    }
    return rsl/(N*p.sigmaT*dx);
}

//calcul de la solution analytique sans scattering, cas homogène 
//et source ponctuelle
double density_no_scattering_homog_point(point p){
    if (p.get_mu()==0) {if (p.get_x()==0){return 1;} return 0;}
    return (1/p.get_mu())*exp(-p.sigmaT*p.get_x()/p.get_mu());}



/* CAS HOMOGENE SOURCE UNIFORME */ 
//tirage de départ de neutrons dans un matériau homogène, source uniforme 
//arrivés en (x,mu).
vector_points no_scattering_homog_unif_MC(int N, point p){
    vector_points selection(N);
    for (int i=0; i<N; i++){
        selection.points[i] = deplacement_x(point(p.get_x(), -p.get_mu()));
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source isotrope uniforme et matériau homogène.
double density_no_scattering_homog_unif_MC(point p, int N){
    vector_points selection(N) ;
    selection = no_scattering_homog_unif_MC(N,p);
    double rsl = 0.;
    for (int i=0;i<N;i++){
        if (selection.points[i].get_x()>0 && selection.points[i].get_x()<1) {
            rsl++;
        }
    }
    //print_vector(selection);
    return rsl/(N*abs(p.get_mu()));
}

//calcul de la solution analytique sans scattering, cas homogène 
//et source uniforme : disjonction cas sur signe de mu 
double density_no_scattering_homog_unif(point p){
    if (p.get_mu()>0) {return 1-exp(-p.sigmaT*p.get_x()/p.get_mu());}
    return 1-exp(p.sigmaT*(1-p.get_x())/p.get_mu());}
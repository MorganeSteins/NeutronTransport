#include "MC.hpp"
#include <math.h>
#include<iostream>
#include <vector>
#include <cstdlib>

using namespace std;

//tirage de départ de neutrons arrivés en 1
vector_points no_scattering_MC(int N, double mu){
    vector_points selection(N);
    for (int i=0; i<N; i++){
        selection.points[i] = point(parcourt_x(point(0,mu)),mu); //retour au pt de depart du neutron
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source ponctuelle isotrope en 0
double density_no_scattering_MC(point p, int N){
    vector_points selection ;
    selection = no_scattering_MC(N,p.get_mu());
    double rsl = 0.;
    for (int i=0;i<N;i++){
        if (selection.points[i].get_x()>p.get_x()) {
            rsl++;
        }
    }
    //print_vector(selection);
    return rsl/N;
}

//calcul de la solution analytique sans absorption
double density_no_scattering(point p){
    if (p.get_mu()==0) {return 1;}
    return (1/p.get_mu())*exp(-p.sigmaT*p.get_x()/p.get_mu());}
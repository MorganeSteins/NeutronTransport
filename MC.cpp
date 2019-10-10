#include "MC.hpp"
#include <math.h>
#include<iostream>
#include <vector>
#include <cstdlib>

using namespace std;

/* CAS ABSORBANT HOMOGENE SOURCE PONCTUELLE */
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
    double dx=sqrt(1./N);
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



/* CAS ABSORBANT HOMOGENE SOURCE UNIFORME */ 
//tirage de départ de neutrons dans un matériau homogène, source uniforme 
//arrivés en (x,mu).
vector_points no_scattering_homog_unif_MC(int N, point p){
    vector_points selection(N);
    for (int i=0; i<N; i++){
        selection.points[i] = deplacement_x(point(new_x(), p.get_mu()));
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source isotrope uniforme et matériau homogène.
double density_no_scattering_homog_unif_MC(point p, int N){
    vector_points selection(N) ;
    selection = no_scattering_homog_unif_MC(N,p);
    double dx = sqrt(1./N);
    double rsl = 0.;
    for (int i=0;i<N;i++){
        if (selection.points[i].get_x()>p.get_x() && selection.points[i].get_x()<p.get_x()+dx) {
            rsl++;
        }
    }
    //print_vector(selection);
    return rsl/(N*p.sigmaT*dx);
}

//calcul de la solution analytique sans scattering, cas homogène 
//et source uniforme : disjonction cas sur signe de mu 
double density_no_scattering_homog_unif(point p){
    if (p.get_mu()==0) {return 1;}
    if (p.get_mu()<0) {return (1./ p.sigmaT)*(1-exp(-p.sigmaT*(p.get_x()-1)/p.get_mu()));}
    return (1./ p.sigmaT)*(1-exp(-p.sigmaT*p.get_x()/p.get_mu()));}



/* CAS DIFFUSANT HOMOGENE SOURCE UNIFORME */
/* Génération des nouveaux points (x_n+1,mu_n) à partir de (x_n,mu_n) 
 Ne conserve que ceux qui restent dans l'intervalle [0,1] */
vector_points move_scattering(vector_points start){
    vector_points final(start.points.size());
    double new_x_;
    int stop=0;
    int count=0;
    for (int i=0;i<start.points.size();i++){
        new_x_ = deplacement_x(start.points[i]).get_x();
        stop = do_I_stop(point::sigmaS,point::sigmaA);
        if (new_x_>=0 && new_x_<=1 && stop==1){
            final.points[i] = point(new_x_,new_mu());
            count ++;
        }
    }
    final.points.resize(count);
    return final;
}


double density_tilda_scattering_homg_unif_MC(int N, double x, int max_iter, double epsilon ){
    vector_points selection(N);
    double dx = sqrt(1./N);
    double phi=0.; // valeur de l'estimateur
    double freq = 0.;
    double phi_old = 1.; //ancien phi
    int step = 0;
    while (abs(phi-phi_old)>epsilon && step<max_iter && selection.points.size()>0){
        phi_old = phi;
        selection = move_scattering(selection); // resize a eu lieu + sortie et absorbés
        for (int i=0;i<selection.points.size();i++) {
            if (selection.points[i].get_x()>x && selection.points[i].get_x()<x+dx){
                freq ++;
            }
        }
        phi += (1./(point::sigmaT*dx))*(freq/(float)N);
        freq = 0.;
        step++;
    }
    return phi;
}

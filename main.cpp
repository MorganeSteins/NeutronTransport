#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include "points.hpp"
#include "aleat_tools.hpp"
#include "MC.hpp"

using namespace std;

int main() {
    srand (static_cast <unsigned> (time(0)));
    
    // Point où on travaille
    double x = 1;
    double mu = 1;
    point p(x,mu);

    // // QUESTION 4
    // // Calcul de la densité au point (x,mu)
    int N = 1e5;
    // cout<<"Calcul de phi q4 au point ("<<x<<","<<mu<<") avec "<<N<<" tirages"<<endl;
    // cout<<"densité par MC "<<density_no_scattering_homog_point_MC(point(x,mu), N)<<endl;
    // cout<<"densité réelle "<<density_no_scattering_homog_point(point(x,mu))<<endl;
    
    // //calcul de la densité en nb_points intervalles et sauvée dans un fichier txt pour les plots
    int nb_points = 25;
    // N = 1e4;
    // cout<<"N = "<<N<<" et on étudie "<<nb_points<<" intervalles."<<endl<<endl;
    // density_segment_no_scattering_homog_point_MC(N,nb_points,1);


    // // QUESTION 5
    // // Calcul de la densité au point (x,mu)
    // N = 1e6;
    // cout<<"Calcul de phi q5 au point ("<<x<<","<<mu<<") avec "<<N<<" tirages"<<endl;
    // cout<<"densité par MC "<<density_no_scattering_homog_unif_MC(point(x,mu), N)<<endl;
    // cout<<"densité réelle "<<density_no_scattering_homog_unif(point(x,mu))<<endl;

    // //calcul de la densité en nb_points intervalles et sauvée dans un fichier txt pour les plots
    // nb_points = 25;
    // N = 1e5;
    // cout<<"N = "<<N<<" et on étudie "<<nb_points<<" intervalles."<<endl<<endl;
    // density_segment_no_scattering_homog_unif_MC(N,nb_points,p);
    


    // QUESTION 8
    //Cette fois ci on ne calcule que le flux moyenné en mu
    nb_points = 25;
    N = 1e5;
    int max_iter = 1e5;
    double epsilon = 1e-10;

    //density_scattering_homog_unif_MC(1e7,p,max_iter,epsilon);
    density_tilda_segment_scattering_homg_unif_MC(N,nb_points, max_iter, epsilon);



    return 0;
}

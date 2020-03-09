#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "points.hpp"
#include "aleat_tools.hpp"
// #include "MC.hpp"
#include "Mc_iter.hpp"
// #include "matplotlibcpp.h"
// #include "deterministe.hpp"
#include "cases.hpp"

#include <chrono>
#include <unistd.h>

using namespace std;
// namespace plt = matplotlibcpp;

int main()
{
    srand(static_cast<unsigned>(time(0)));

    // Point où on travaille
    double x = 1;
    double mu = 1;
    point p(x, mu);

    // // QUESTION 4
    // // Calcul de la densité au point (x,mu)
    int N = 1e7;
    // int N= 100 ;
    int max_iter = 100;
    // cout<<"Calcul de phi q4 au point ("<<x<<","<<mu<<") avec "<<N<<" tirages"<<endl;
    // cout<<"densité par MC "<<density_no_scattering_homog_point_MC(point(x,mu), N)<<endl;
    // cout<<"densité réelle "<<density_no_scattering_homog_point(point(x,mu))<<endl;

    // //calcul de la densité en nb_points intervalles et sauvée dans un fichier txt pour les plots
    int nb_points = 25;

    cout << "\n ----TEST poly---- \n";
    double data1 = density_no_scattering_MC<source_uniforme> (mu, x , N);
    double data2 = density_no_scattering_homog_unif_MC_iter( mu,x, N) ;
    cout<<"Résultat "<<data1<<" et autre "<<data2<<"\n" ;

    // N = 1e4;
    // cout<<"N = "<<N<<" et on étudie "<<nb_points<<" intervalles."<<endl<<endl;
    // cout << " ----TEST 1---- \n";
    // auto start = chrono::steady_clock::now();
    // auto data1 = density_no_scattering_homog_point_MC(p, N);
    // auto end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    // start = chrono::steady_clock::now();
    // auto data2 = density_no_scattering_homog_point_MC_iter(p, N);
    // end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";
    // cout << "Error comparison " << data2 - data1 << "\n\n";

    // cout << "\n ----TEST 2---- \n";
    // start = chrono::steady_clock::now();
    // data1 = density_tilda_no_scattering_homog_point_MC(p, N);
    // end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    // start = chrono::steady_clock::now();
    // data2 = density_tilda_no_scattering_homog_point_MC_iter(p, N);
    // end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    // start = chrono::steady_clock::now();
    // auto data3 = density_tilda_no_scattering_homog_point_MC_iter(x, N);
    // end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n\n";
    // cout << "Error comparison " << data2 - data1 << "\n";
    // cout << "Error comparison " << data3 - data1 << "\n\n";

    // cout << "\n ----TEST 3---- \n";
    // start = chrono::steady_clock::now();
    // auto liste1 = density_segment_no_scattering_homog_unif_MC(N, nb_points, p);
    // end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    // start = chrono::steady_clock::now();
    // auto liste2 = density_segment_no_scattering_homog_unif_MC_iter(N, nb_points, x, mu);
    // end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";
    // cout << "Error comparison " << liste2[0] - liste1[0] << "\n\n";

    // cout << "\n ----TEST 4---- \n";
    // auto start = chrono::steady_clock::now();
    // auto liste1 = density_tilda_segment_scattering_homg_unif_MC( N,  nb_points,  max_iter) ;
    // auto end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    // start = chrono::steady_clock::now();
    // vector<double> liste2 = density_tilda_segment_scattering_homg_unif_MC_iter( N,  nb_points,  max_iter) ;
    // end = chrono::steady_clock::now();
    // cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";
    // cout << "Error comparison " << liste2[0] - liste1[0] << "\n\n";

    
    
    // vector<double> x_plot(nb_points) ;
    // for (int i=0;i<N;i++) x_plot[i] = i ;

    // plt::plot(liste1) ;
    // plt::plot(liste2) ;
    // plt::legend() ; 
    // plt::show();

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

    //density_scattering_homog_unif_MC(1e7,p,max_iter,epsilon);
    // density_tilda_segment_scattering_homg_unif_MC(N,nb_points, max_iter);
    // density_tilda_segment_scattering_woodcock_unif_MC(N, nb_points, max_iter, epsilon);

    return 0;
}

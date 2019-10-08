#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <time.h>
#include "points.hpp"
#include "aleat_tools.hpp"
#include "MC.hpp"

using namespace std;

int main() {
    srand (static_cast <unsigned> (time(0)));
    
    //print_vector(sol);
    double x = 0.1;
    double mu = 0.5;
    double integrale = 0.;
    cout<<"On travaille à x="<<x<<" et mu="<<mu<<endl;

    // test question 4
    double density = density_no_scattering_homog_point_MC(point(x,mu), 1000000);
    cout<<"densité par MC "<<density<<endl;
    cout<<"densité réelle "<<density_no_scattering_homog_point(point(x,mu))<<endl<<endl;
    
    // test question 5
    double density_2 = density_no_scattering_homog_unif_MC(point(x,mu), 1000000);
    cout<<"densité par MC "<<density_2<<endl;
    cout<<"densité réelle "<<density_no_scattering_homog_unif(point(x,mu))<<endl<<endl;
    
    //for (int i=0;i<1000;i++) {integrale += (15./1000.)*density_no_scattering_homog_unif(point((float) i*15/1000,mu));}
    //cout<<integrale<<endl;
    //for (int i=0;i<11;i++) {cout<<(float)i/10<<" vaut "<<density_no_scattering_homog_unif_MC(point(float(i)/10,mu),10000)<<endl;}
    /* test du flux non collisionné : intégrale somme à 1*/

    // somme des valeurs pour approcher l'intégrale MC : doit faire 1. De -7 à 7.
    for (int i=0;i<1000;i++) {integrale += (15./1000.)*density_no_scattering_homog_unif_MC(point((float) i*15/1000-7,mu),10000);}
    cout<<integrale<<endl;
    return 0;
}

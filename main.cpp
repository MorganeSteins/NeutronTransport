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
    
    vector_points sol;
    sol = no_scattering_MC(2,0.);
    //print_vector(sol);
    double x = 1;
    double mu = 0.25;
    double integrale = 0.;
    cout<<"On travaille à x="<<x<<" et mu="<<mu<<endl;


    /* test du flux non collisionné : intégrale somme à 1*/
    /*double density = density_no_scattering_MC(point(x,mu), 100000);
    cout<<"densité par MC "<<density<<endl;
    cout<<"densité réelle "<<density_no_scattering(point(x,mu))<<endl<<endl;
    for (int i=0;i<1000;i++) {integrale += (15./1000.)*density_no_scattering_MC(point((float) i*15/1000,mu),10000);}
    cout<<integrale<<endl;*/
    return 0;
}

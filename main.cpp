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

    //double density = density_no_scattering_MC(point(0.5,0.5), 10000);
    //cout<<"densité par MC "<<density<<endl;
    cout<<"densité réelle mu=0.1  "<<density_no_scattering(point(0.5,0.1))<<endl;
    cout<<"densité réelle mu=0.5  "<<density_no_scattering(point(0.5,0.5))<<endl;
    cout<<"densité réelle mu=1  "<<density_no_scattering(point(0.5,1))<<endl;
    for (int i=0;i<10;i++){cout<<(float) i/10<<" donne "<<density_no_scattering(point(1,(float) i/10))<<endl;}
    return 0;
}

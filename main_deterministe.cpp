#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include "points.hpp"
#include "aleat_tools.hpp"
#include "MC.hpp"
#include "deterministe.hpp"

int main(int argc, char *argv[]) {
    double x,mu;
    if (argc != 3)
    {
        printf("You need to input 2 variables:x and mu, we are using default value 0,1 \n");
        x  = 0;
        mu = 1.;
    }

    else {
        x  = atof(argv[1]);
        mu = atof(argv[2]);
    }
    cout<<argv[1]<<"  "<<argv[2]<<endl;
    int Nx=5;
    int Nmu = 2;
    double dx = 1./Nx;
    vector<double> Q(Nx);
    vector<double> sigmaT(Nx);
    for (int i=0;i<Nx;i++){
        Q[i] = 0;
        if (i*dx>=0.3 && i*dx<0.7) {sigmaT[i]=3;}
        else {sigmaT[i]=1;}
    }

    vector<double> phi = IS(Nx,Nmu,0.01,100, Q,sigmaT);
    // for (int i=0;i<phi.size();i++) {cout<<"x="<<i*dx<<" "<<sigmaT[i]<<endl;}
    return 0;
}
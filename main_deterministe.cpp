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
    int Nx,Nmu;
    if (argc != 3)
    {
        printf("You need to input 2 variables:Nx and Nmu (EVEN NUMBER!!), we are using default value 5 and 2 \n");
        Nx  = 5;
        Nmu = 2;
    }

    else {
        Nx  = atof(argv[1]);
        Nmu = atof(argv[2]);
    }

    double dx = 1./Nx;
    vector<double> Q(Nx);
    vector<double> sigmaT(Nx);
    for (int i=0;i<Nx;i++){
        Q[i] = 1;
        sigmaT[i] = 1;
        // if (i*dx>=0.3 && i*dx<0.7) {sigmaT[i]=3;}
        // else {sigmaT[i]=1;}
    }

    vector<double> phi = IS(Nx,Nmu,0.01,5, Q,sigmaT);
    // for (int i=0;i<phi.size();i++) {cout<<"x="<<i*dx<<" "<<sigmaT[i]<<endl;}
    return 0;
}
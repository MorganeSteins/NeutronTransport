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

    double epsilon = 0.001;
    double dx = 1./Nx;
    vector<double> Q(Nx);
    vector<double> sigmaT(Nx);
    vector<double> sigmaS(Nx);
    double sa=1;
    for (int i=0;i<Nx;i++){
        Q[i] = 1*epsilon;
        sigmaT[i] = 1/epsilon;
        sigmaS[i] = sigmaT[i]-sa*epsilon;
        // if (i*dx>=0.3 && i*dx<0.7) {sigmaT[i]=3;}
        // else {sigmaT[i]=1;}
    }

    double mu = 1.;
    double nu=1e-7;
    // vector<double> phi;
    // vector<double> phi = IS_iteration(Nx,mu,1./mu,Q,sigmaT);
    vector<double> phi = IS(Nx,Nmu,nu,500000, Q,sigmaT,sigmaS);
    ofstream fichier("Data/phi_q13_"+to_string(Nx)+"_"+to_string(Nmu)+"_epsilon0001_sa.txt", ios::out | ios::trunc);
    for (int i=0;i<phi.size();i++){
        fichier<<phi[i]<<",";
    }
    fichier.close();
    cout<<"Saved in Data/phi_q13_"+to_string(Nx)+"_"+to_string(Nmu)+"_epsilon01_sa.txt"<<endl;
    // for (int i=0;i<phi.size();i++) {cout<<"x="<<i*dx+dx/2<<" "<<phi[i]<<endl;}
    // for (int i=10;i<=1000000;i*=10){
    //     cout<<"Nx="<<i<<endl;
    //     Nx=i/2;
    //     dx = 1./Nx;
    //     vector<double> Q(Nx);
    //     vector<double> sigmaT(Nx);
    //     for (int i=0;i<Nx;i++){
    //         Q[i] = 0;
    //         sigmaT[i] = 1;
    //     }   
    //     phi = IS_iteration(i/2,mu,1./mu,Q,sigmaT);
    //     ofstream fichier("Data/phi_q9_"+to_string(i/2)+"_"+to_string(abs((mu)))+"test.txt", ios::out | ios::trunc);
    //     for (int i=0;i<phi.size();i++){
    //         fichier<<phi[i]<<",";
    //     }
    //     fichier.close();
        
    //     Nx=i;
    //     dx = 1./Nx;
    //     Q.resize(Nx);
    //     sigmaT.resize(Nx);
    //     for (int i=0;i<Nx;i++){
    //         Q[i] = 0;
    //         sigmaT[i] = 1;
    //     }   
    //     phi = IS_iteration(i,mu,1./mu,Q,sigmaT);
    //     IS_iteration(i,mu,1./mu,Q,sigmaT);
    //     ofstream fichier2("Data/phi_q9_"+to_string(i)+"_"+to_string(abs((mu)))+"test.txt", ios::out | ios::trunc);
    //     for (int i=0;i<phi.size();i++){
    //         fichier2<<phi[i]<<",";
    //     }
    //     fichier2.close();
    // }

    return 0;
}
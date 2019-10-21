#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "deterministe.hpp"

vector<double> solveur(int Nx, double mu, double source){
    vector<double> X(Nx+1);
    vector<double> phi(Nx+1);
    double dx = 1./Nx;
    vector<double> Q(Nx);
    int sign_mu = sgn(mu); //signe de mu

    //initialisation de x et Q
    for (int i=0;i<=Nx;i++){
        X[i] = i*dx;
        Q[i] = 0;
    }

    int debut = Nx*(1-sign_mu)/2;
    int fin = Nx*(1+sign_mu)/2;

    phi[debut] = source;
    cout<<debut<<" a "<<fin<<endl;

    for (int i=debut+1;i<=fin;i+=sign_mu){
        phi[i] = (2*dx*Q[i-1]+(2*abs(mu)-dx*point::sigmaT)*phi[i-1])/(2*abs(mu)+dx*point::sigmaT);
    }


    ofstream fichier("Data/phi_q9_"+to_string(Nx)+".txt", ios::out | ios::trunc);
    for (int i=0;i<phi.size();i++){
        fichier<<phi[i]<<",";
    }
    fichier.close();
    return phi;
}
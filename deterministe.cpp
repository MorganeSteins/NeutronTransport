#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "deterministe.hpp"

vector<double> solveur(int Nx, double mu, double source, vector<double> Q,vector<double> sigmaT){
    vector<double> X(Nx+1);
    vector<double> phi(Nx+1);
    double dx = 1./Nx;
    int sign_mu = sgn(mu); //signe de mu

    //initialisation de x et Q
    for (int i=0;i<=Nx;i++){
        X[i] = i*dx;
    }

    int debut = Nx*(1-sign_mu)/2;
    int fin = Nx*(1+sign_mu)/2;

    phi[debut] = abs(source);
    cout<<debut<<" a "<<fin<<endl;

    for (int i=debut;i!=fin;i+=sign_mu){
        phi[i+sign_mu] = (2*dx*Q[i]+(2*abs(mu)-dx*sigmaT[i])*phi[i])/(2*abs(mu)+dx*sigmaT[i]);
        // cout<<"i="<<i<<"  "<<phi[i+sign_mu]<<endl;
    }


    ofstream fichier("Data/phi_q11_"+to_string(Nx)+"_"+to_string(abs((mu)))+".txt", ios::out | ios::trunc);
    for (int i=0;i<phi.size();i++){
        fichier<<phi[i]<<",";
    }
    fichier.close();
    return phi;
}
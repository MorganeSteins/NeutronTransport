#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include "src/deterministe/deterministe.hpp"

int main(int argc, char *argv[])
{
    //Reading inputs
    int Nx, Nmu;
    if (argc != 3)
    {
        printf("You need to input 2 variables:Nx and Nmu (EVEN NUMBER!!), we are using default value 5 and 2 \n");
        Nx = 5;
        Nmu = 2;
    }

    else
    {
        Nx = atof(argv[1]);
        Nmu = atof(argv[2]);
    }

    // Initialize parameters
    double epsilon = 1;
    double dx = 1. / Nx;
    vector<double> Q(Nx);
    vector<double> sigmaT(Nx);
    vector<double> sigmaS(Nx);
    double sa = 0;
    for (int i = 0; i < Nx; i++)
    {
        Q[i] = 1 * epsilon;
        sigmaT[i] = 1 / epsilon;
        sigmaS[i] = sigmaT[i] - sa * epsilon;
    }

    double mu = 1.;
    double eta = 1e-7;
    int max_iter = 1000000;

    // Compute Source iteration
    vector<double> phi = IS(Nx, Nmu, eta, max_iter, Q, sigmaT, sigmaS);

    // Saves output
    string save_filename = "Data/phi_q15_" + to_string(Nx) + "_" + to_string(Nmu) + "_" + to_string(epsilon) + ".txt";
    ofstream fichier(save_filename, ios::out | ios::trunc);
    for (int i = 0; i < phi.size(); i++)
    {
        fichier << phi[i] << ",";
    }
    fichier.close();

    // User message
    cout << save_filename << endl;
    return 0;
}
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include "deterministe.hpp"

#include "eigen/Eigen/IterativeLinearSolvers"
#include "eigen/Eigen/Sparse"

int main(int argc, char *argv[])
{
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

    double epsilon = 0.1, sa = 1;
    double dx = 1. / Nx;
    Eigen::VectorXd Q(Nx), sigmaT(Nx + 1), sigmaS(Nx + 1);
    Q.setConstant(1 * epsilon);
    sigmaT.setConstant(1. / epsilon);
    sigmaS = sigmaT - sa * epsilon * Eigen::VectorXd::Ones(Nx + 1);
    // sigmaS = sigmaT;

    // cout << sigmaT;

    double nu = 1e-3;
    Eigen::VectorXd phi = Fast_IS(Nx, Nmu, 1e-7, 1000, Q, sigmaT, sigmaS);

    ofstream fichier("Data/phi_fast_" + to_string(Nx) + "_" + to_string(Nmu) + "_epsilon" + to_string(epsilon) + ".txt", ios::out | ios::trunc);
    fichier << phi;
    fichier.close();

    cout << "Saved in Data/phi_fast_" + to_string(Nx) + "_" + to_string(Nmu) + "_epsilon" + to_string(epsilon) + ".txt" << endl;

    return 0;
}
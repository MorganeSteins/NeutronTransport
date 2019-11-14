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
    int Nx, mu;
    if (argc != 3)
    {
        printf("You need to input 2 variables:Nx and mu, we are using default value 2 and 1/2 \n");
        Nx = 2;
        mu = 0.5;
    }
    else
    {
        Nx = atof(argv[1]);
        mu = atof(argv[2]);
    }

    double flux_entrant = 0;

    vector<double> Q;
    for (size_t i = 0; i < Nx; i++)
    {
        Q.push_back(1.);
    }
    auto sigmaT = Q;

    auto res1 = IS_iteration(Nx, mu, flux_entrant, Q, sigmaT);

    auto QQ = Eigen::VectorXd::Constant(Nx, 1);
    auto sigmaTT = QQ;

    auto res2 = IS_iteration(Nx, mu, flux_entrant, QQ, sigmaTT);

    cout << "res 1 = ";
    for (size_t i = 0; i < Nx; i++)
        cout << res1[i] << " ";
    cout << endl;

    cout << "res 2 = " << res2.transpose() << endl;

    return 0;
}
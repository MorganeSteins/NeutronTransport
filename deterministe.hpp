#ifndef DETERMINISTE_INCLUDED
#define DETERMINISTE_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include "points.hpp"

// eigen inclusion
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/IterativeLinearSolvers"

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

vector<double> IS_iteration(int Nx, double mu, double flux_entrant, vector<double> Q,vector<double> sigmaT);
vector<double> IS(int Nx, int Nmu, double epsilon, int iter_max, vector<double> S, vector<double> sigmaT,vector<double> sigmaS);

Eigen::VectorXd IS_iteration(int Nx, double mu, double flux_entrant, Eigen::VectorXd Q, Eigen::VectorXd sigmaT);
Eigen::VectorXd IS_iteration_phi(int Nx, double mu, double flux_entrant, Eigen::VectorXd Q, Eigen::VectorXd sigmaT);

Eigen::VectorXd Fast_IS(int Nx, int Nmu, double epsilon, int iter_max, Eigen::VectorXd S, Eigen::VectorXd sigmaT, Eigen::VectorXd sigmaS);

#endif // DETERMINISTE_INCLUDED
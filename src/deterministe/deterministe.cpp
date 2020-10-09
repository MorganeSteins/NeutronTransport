#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "deterministe.hpp"

#include "../../eigen/Eigen/IterativeLinearSolvers"
#include "../../eigen/Eigen/Sparse"

vector<double> IS_iteration(int Nx, double mu, double flux_entrant, vector<double> Q, vector<double> sigmaT)
{
    vector<double> phi_t(Nx);
    double dx = 1. / Nx;
    double phi_plus;
    int sign_mu = sgn(mu); //signe de mu

    double phi_moins = std::abs(flux_entrant);

    if (sign_mu > 0)
    {
        for (int i = 0; i < Nx; i++)
        {
            phi_plus = (2 * dx * Q[i] + (2 * std::abs(mu) - dx * sigmaT[i]) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT[i]);
            phi_t[i] = (1. / 2) * (phi_plus + phi_moins);
            phi_moins = phi_plus;
        }
    }
    else
    {
        for (int i = Nx - 1; i >= 0; i--)
        {
            phi_plus = (2 * dx * Q[i] + (2 * std::abs(mu) - dx * sigmaT[i]) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT[i]);
            phi_t[i] = (1. / 2) * (phi_plus + phi_moins);
            phi_moins = phi_plus;
        }
    }
    return phi_t;
}

//solveur IS
vector<double> IS(int Nx, int Nmu, double epsilon, int iter_max, vector<double> S, vector<double> sigmaT, vector<double> sigmaS)
{
    vector<double> Q(Nx); //source à initialiser selon le pb
    vector<double> X(Nx + 1);
    vector<double> MU(Nmu);
    vector<double> phi(Nx);
    vector<double> Q2(Nx);

    //initialisation de X
    double dx = 1. / Nx;
    for (int i = 0; i <= Nx; i++)
        X[i] = i * dx;

    //initialisation de MU
    double dmu = 2. / Nmu;
    double w = dmu; //poids de quadrature constant
    for (int i = 0; i < Nmu; i++)
        MU[i] = i * dmu - 1 + dmu / 2;

    //initialisation de Q
    for (int i = 0; i < Nx; i++)
    {
        Q[i] = S[i];
        Q2[i] = 0;
    }
    double err = 1.; // erreur que l'on initialise grande
    double norm_Q = 1.;
    int n_iter = 0; //compteur itération
    while (err / norm_Q > epsilon * epsilon && n_iter < iter_max)
    {
        for (int n = 0; n < Nmu; n++)
        {
            phi = IS_iteration(Nx, MU[n], 0., Q, sigmaT); // ATTENTION FLUX ENTRANT
            for (int i = 0; i < Nx; i++)
            {
                Q2[i] += (1. / 2) * sigmaS[i] * w * phi[i];
            }
        }
        err = 0.;
        norm_Q = 0.;

        for (int i = 0; i < Nx; i++)
        {
            auto tmp = (Q[i] - (Q2[i] + S[i]));
            err += tmp * tmp;
            Q[i] = Q2[i] + S[i];
            norm_Q += Q[i] * Q[i];
            Q2[i] = 0;
        }
        n_iter++;
    }
    cout << "Convergé en " << n_iter << " itérations avec une erreur " << sqrt(err / norm_Q) << endl;
    cout << endl;

    vector<double> phi_final(Nx, 0.);
    for (int n = 0; n < Nmu; n++)
    {
        phi = IS_iteration(Nx, MU[n], 0., Q, sigmaT); // ATTENTION FLUX ENTRANT NUL
        for (int i = 0; i < Nx; i++)
            phi_final[i] += (1. / 2) * w * phi[i];
    }
    return phi_final;
}

Eigen::VectorXd IS_iteration(int Nx, double mu, double flux_entrant, Eigen::VectorXd Q, Eigen::VectorXd sigmaT)
{
    using namespace Eigen;

    VectorXd phi_t(Nx);
    double dx = 1. / Nx;
    double phi_plus;
    int sign_mu = sgn(mu);

    double phi_moins = std::abs(flux_entrant);

    if (sign_mu > 0)
    {
        for (int i = 0; i < Nx; i++)
        {
            phi_plus = (2 * dx * Q(i) + (2 * std::abs(mu) - dx * sigmaT(i)) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT(i));
            phi_t(i) = (1. / 2.) * (phi_plus + phi_moins);
            phi_moins = phi_plus;
        }
    }
    else
    {
        for (int i = Nx - 1; i >= 0; i--)
        {
            phi_plus = (2 * dx * Q(i) + (2 * std::abs(mu) - dx * sigmaT(i)) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT(i));
            phi_t(i) = (1. / 2.) * (phi_plus + phi_moins);
            phi_moins = phi_plus;
        }
    }

    return phi_t;
}

Eigen::VectorXd IS_iteration_phi(int Nx, double mu, double flux_entrant, Eigen::VectorXd Q, Eigen::VectorXd sigmaT)
{
    using namespace Eigen;

    VectorXd phi(Nx + 1);
    double dx = 1. / Nx;
    int sign_mu = sgn(mu);

    if (sign_mu > 0)
    {
        double phi_moins = phi(0) = std::abs(flux_entrant);
        for (size_t i = 1; i < Nx + 1; i++)
        {
            phi(i) = phi_moins = (2 * dx * Q(i) + (2 * std::abs(mu) - dx * sigmaT(i)) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT(i));
        }
    }
    else
    {
        double phi_moins = phi(Nx) = std::abs(flux_entrant);
        for (size_t i = Nx - 1; i >= 0; i--)
        {
            phi(i) = phi_moins = (2 * dx * Q(i) + (2 * std::abs(mu) - dx * sigmaT(i)) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT(i));
        }
    }

    return phi;
}

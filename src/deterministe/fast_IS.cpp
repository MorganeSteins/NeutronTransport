#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "deterministe.hpp"

#include "../../eigen/Eigen/IterativeLinearSolvers"
#include "../../eigen/Eigen/Sparse"

Eigen::VectorXd Fast_IS(int Nx, int Nmu, double epsilon, int iter_max, const Eigen::VectorXd S, const Eigen::VectorXd sigmaT, const Eigen::VectorXd sigmaS)
{
    using namespace Eigen;

    VectorXd Q(Nx), Q2(Nx); // source à initialiser selon le pb

    VectorXd X(Nx + 1);          // subdivision spatiale
    VectorXd MU(Nmu);            // subdivision angulaire
    VectorXd phi_tilde(Nx);      // solution courante
    VectorXd phi_demi_tilde(Nx); // solution au temps t+1/2
    VectorXd q_demi(Nx);         //
    VectorXd L(Nx - 1);          // Second membre du système DSA
    VectorXd F(Nx + 1);          // solution du système DSA

    //initialisation de X
    double dx = 1. / Nx;
    for (int i = 0; i <= Nx; i++)
        X(i) = i * dx;

    //initialisation de MU
    double dmu = 2. / Nmu;
    double w = dmu; //poids de quadrature constant
    for (int i = 0; i < Nmu; i++)
        MU(i) = (i + 1. / 2.) * dmu - 1;

    // sigma_s tilde and sigma_t tilde
    VectorXd sigmaS_tilde = 0.5 * (sigmaS.head(Nx) + sigmaS.tail(Nx));
    VectorXd sigmaT_tilde = 0.5 * (sigmaT.head(Nx) + sigmaT.tail(Nx));

    // Définition de la matrice de diffusion
    SparseMatrix<double> A(Nx - 1, Nx - 1);
    for (size_t i = 0; i < Nx - 1; i++)
    {
        A.coeffRef(i, i) = 2. / (3. * dx * sigmaT(i)) + (sigmaT(i) - sigmaS(i)) * dx * 2. / 3.;
        if (i < Nx - 2)
        {
            A.coeffRef(i + 1, i) = -1. / (3 * dx * sigmaT(i)) + (sigmaT(i) - sigmaS(i)) * dx / 6.;
            A.coeffRef(i, i + 1) = A.coeffRef(i + 1, i);
        }
    }

    // initialisation du solveur
    SimplicialLDLT<SparseMatrix<double>> solver_ldlt;
    solver_ldlt.compute(A);

    // matrice q->L
    SparseMatrix<double> q2L(Nx + 1, Nx);
    for (size_t i = 0; i < Nx; i++)
    {
        q2L.coeffRef(i, i) = 1;
        q2L.coeffRef(i + 1, i) = 1;
    }
    q2L *= dx / 2.;

    //initialisation de Q
    Q = S;
    Q2.setZero(Nx);

    double err = 1.;    // erreur que l'on initialise grande pour entre dans la boucle
    double norm_Q = 1.; //
    int n_iter = 0;     //compteur itération

    while (err / norm_Q > epsilon * epsilon && n_iter < iter_max)
    {
        // solveur diamant
        phi_demi_tilde.setZero(Nx);
        for (size_t n = 0; n < Nmu; n++)
        {
            phi_demi_tilde += w / 2 * IS_iteration(Nx, MU(n), 0., Q, sigmaT_tilde);
        }

        // DSA
        q_demi = (phi_demi_tilde - phi_tilde).cwiseProduct(sigmaS_tilde);

        L = (q2L * q_demi).segment(1, Nx - 1);
        F.segment(1, Nx - 1) = solver_ldlt.solve(L);

        F(0) = 0.;
        F(Nx) = 0.;

        // nouvelle solution
        phi_tilde = phi_demi_tilde + 0.5 * (F.head(Nx) + F.tail(Nx));

        // nouvelle source
        Q2 = phi_tilde.cwiseProduct(sigmaS_tilde);
        err = (Q - Q2 - S).dot(Q - Q2 - S);
        norm_Q = Q.dot(Q);
        Q = Q2 + S;
        Q2.setZero(Nx);

        n_iter++;
    }
    err = sqrt(err);
    cout << "Converged in " << n_iter << " iterations with error " << err << endl;
    cout << endl;

    return phi_tilde;
}
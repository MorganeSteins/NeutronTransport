#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "deterministe.hpp"

#include "eigen/Eigen/IterativeLinearSolvers"
#include "eigen/Eigen/Sparse"

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
            // cout<<"phi moins = "<<phi_moins<<" et phi_plus = "<<phi_plus<<endl;
            phi_moins = phi_plus;
            // cout << "i=" << i << "  " << phi_t[i] << endl;
        }
    }
    else
    {
        for (int i = Nx - 1; i >= 0; i--)
        {
            phi_plus = (2 * dx * Q[i] + (2 * std::abs(mu) - dx * sigmaT[i]) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT[i]);
            phi_t[i] = (1. / 2) * (phi_plus + phi_moins);
            // cout<<"phi moins = "<<phi_moins<<" et phi_plus = "<<phi_plus<<endl;
            phi_moins = phi_plus;
            // cout << "i=" << i << "  " << phi_t[i] << endl;
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
        // cout<<" difference entre 2 itérations "<<sqrt(err)<<endl;
        n_iter++;
    }
    cout << "Convergé en " << n_iter << " itérations avec une erreur " << sqrt(err) << endl;
    cout << endl;

    vector<double> phi_final(Nx, 0.);
    for (int n = 0; n < Nmu; n++)
    {
        phi = IS_iteration(Nx, MU[n], 0., Q, sigmaT); // ATTENTION FLUX ENTRANT
        for (int i = 0; i < Nx; i++)
            phi_final[i] += (1. / 2) * w * phi[i];
    }

    return phi_final;
}

// vector<double> IS_iteration_phi(int Nx, double mu, double flux_entrant, vector<double> Q,vector<double> sigmaT){
//     vector<double> X(Nx+1);
//     vector<double> phi(Nx+1);
//     double dx = 1./Nx;
//     int sign_mu = sgn(mu); //signe de mu

//     //initialisation de x et Q
//     for (int i=0;i<=Nx;i++){
//         X[i] = i*dx;
//     }

//     int debut = Nx*(1-sign_mu)/2;
//     int fin = Nx*(1+sign_mu)/2;

//     phi[debut] = abs(flux_entrant);
//     cout<<debut<<" a "<<fin<<endl;

//     for (int i=debut;i!=fin;i+=sign_mu){
//         phi[i+sign_mu] = (2*dx*Q[i]+(2*abs(mu)-dx*sigmaT[i])*phi[i])/(2*abs(mu)+dx*sigmaT[i]);
//         // cout<<"i="<<i<<"  "<<phi[i+sign_mu]<<endl;
//     }

//     ofstream fichier("Data/phi_q11_"+to_string(Nx)+"_"+to_string(abs((mu)))+".txt", ios::out | ios::trunc);
//     for (int i=0;i<phi.size();i++){
//         fichier<<phi[i]<<",";
//     }
//     fichier.close();
//     return phi;
// }

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
            // cout<<"phi moins = "<<phi_moins<<" et phi_plus = "<<phi_plus<<endl;
            phi_moins = phi_plus;
            // cout << "i=" << i << "  " << phi_t[i] << endl;
        }
    }
    else
    {
        for (int i = Nx - 1; i >= 0; i--)
        {
            phi_plus = (2 * dx * Q(i) + (2 * std::abs(mu) - dx * sigmaT(i)) * phi_moins) / (2 * std::abs(mu) + dx * sigmaT(i));
            phi_t(i) = (1. / 2.) * (phi_plus + phi_moins);
            // cout<<"phi moins = "<<phi_moins<<" et phi_plus = "<<phi_plus<<endl;
            phi_moins = phi_plus;
            // cout << "i=" << i << "  " << phi_t[i] << endl;
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

Eigen::VectorXd Fast_IS(int Nx, int Nmu, double epsilon, int iter_max, const Eigen::VectorXd S, const Eigen::VectorXd sigmaT, const Eigen::VectorXd sigmaS)
{
    using namespace Eigen;

    VectorXd Q(Nx), Q2(Nx); // source à initialiser selon le pb

    VectorXd X(Nx + 1);          // subdivision spatiale
    VectorXd MU(Nmu);            // subdivision angulaire
    VectorXd phi_tilde(Nx);      // solution courante
    VectorXd phi_demi_tilde(Nx); // solution au temps t+1/2
    VectorXd q_demi(Nx);         //
    VectorXd L(Nx + 1);          // Second membre du système DSA
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
    // A.coeffRef(0, 0) = (sigmaT(0) - sigmaS(0)) * dx * 1. / 3. - 1. / (3. * dx * sigmaT(0));
    // A.coeffRef(Nx - 2, Nx - 2) = (sigmaT(Nx) - sigmaS(Nx)) * dx * 1. / 3. - 1. / (3. * dx * sigmaT(Nx));

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

        cout << "\nphi_demi_tilde = \n";
        cout << phi_demi_tilde << endl;

        // DSA
        q_demi = (phi_demi_tilde - phi_tilde).cwiseProduct(sigmaS_tilde);

        cout << "\nq_demi = \n";
        cout << q_demi << endl;

        L = (q2L * q_demi).segment(1, Nx - 1);
        F.segment(1, Nx - 1) = solver_ldlt.solve(L);

        cout << "\nL = \n";
        cout << L << endl;

        F(0) = 0.;
        F(Nx) = 0.;

        cout << "\nF = \n";
        cout << F << endl;

        // nouvelle solution
        phi_tilde = phi_demi_tilde + 0.5 * (F.head(Nx) + F.tail(Nx));

        cout << "\nphi_tilde = \n";
        cout << phi_tilde << endl;

        // nouvelle source
        Q2 = phi_tilde.cwiseProduct(sigmaS_tilde);
        err = (Q - Q2 - S).dot(Q - Q2 - S);
        norm_Q = Q.dot(Q);
        Q = Q2 + S;
        Q2.setZero(Nx);

        n_iter++;
    }
    err = sqrt(err);
    cout << "Convergé en " << n_iter << " itérations avec une erreur " << err << endl;
    cout << endl;

    return phi_tilde;
}
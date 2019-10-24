#include "MC.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

/* CAS ABSORBANT HOMOGENE SOURCE PONCTUELLE */
//tirage de départ de neutrons dans un matériau homogène, source ponctuelle en 0.
vector_points no_scattering_homog_point_MC(int N, double mu)
{
    vector_points selection(N);
    for (int i = 0; i < N; i++)
    {
        selection.points[i] = point(parcours_x(point(0, mu)), mu);
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source ponctuelle isotrope en 0 et matériau homogène.
double density_no_scattering_homog_point_MC(point p, int N)
{
    double dx = 1. / 100;
    vector_points selection(N);
    selection = no_scattering_homog_point_MC(N, p.get_mu());
    double rsl = 0.;
    for (int i = 0; i < N; i++)
    {
        if (selection.points[i].get_x() > p.get_x() - dx / 2 && selection.points[i].get_x() < dx / 2 + p.get_x())
        {
            rsl++;
        }
    }
    return rsl / (N * p.sigmaT * dx);
}

// SOURCE PONCTUELLE. Tirage de la densité neutronique pour mu fixé sur un maillage de nb_points de l'intervalle [0,1]. Sauve phi et la fréquence dans des fichiers .txt
vector<double> density_segment_no_scattering_homog_point_MC(int N, int nb_points, double mu)
{

    double dx = 1. / nb_points;                      //largeur des intervalles
    vector_points selection(N);                      //stocke les N tirages aléatoires
    vector<double> freq(nb_points, 0.);              //stocke la fréqence pour chaque intervalle
    selection = no_scattering_homog_point_MC(N, mu); // libre parcours depuis (0,mu)
    int indice = 0;                                  //initialisation

    for (int i = 0; i < selection.points.size(); i++)
    {
        indice = floor(selection.points[i].get_x() / dx); //indice de l'intervalle où est le neutron
        if (indice < nb_points)
        {
            freq[indice] += (1. / (point::sigmaT * dx)) * (1. / (float)N);
        }
    }

    // enregistrement des valeurs de phi
    ofstream fichier("Data/phi_q4_" + to_string(N) + ".txt", ios::out | ios::trunc);
    for (int i = 0; i < freq.size(); i++)
    {
        fichier << freq[i] << ",";
    }
    fichier.close();

    // enregistrement des valeurs de fréquence pour les intervalles de confiance
    ofstream fichier2("Data/freq_q4_" + to_string(N) + ".txt", ios::out | ios::trunc);
    for (int i = 0; i < freq.size(); i++)
    {
        fichier2 << freq[i] * N * point::sigmaT * dx << ",";
    }
    fichier2.close();
    cout << "Data/phi_q4_" + to_string(N) + ".txt" << endl;
    return freq;
}

// calcul de la densité neutronique au point x moyennée en mu avec N tirages
// avec une source ponctuelle isotrope en 0 et matériau homogène.
double density_tilda_no_scattering_homog_point_MC(point p, int N)
{
    double dx = sqrt(1. / N);
    vector_points selection(N);
    selection = no_scattering_homog_point_MC(N, new_mu());
    double rsl = 0.;
    for (int i = 0; i < N; i++)
    {
        if (selection.points[i].get_x() > p.get_x() && selection.points[i].get_x() < dx + p.get_x())
        {
            rsl++;
        }
    }
    return rsl / (2 * N * p.sigmaT * dx);
}

//calcul de la solution analytique sans scattering, cas homogène
//et source ponctuelle
double density_no_scattering_homog_point(point p)
{
    if (p.get_mu() == 0)
    {
        if (p.get_x() == 0)
        {
            return 1;
        }
        return 0;
    }
    return (1 / p.get_mu()) * exp(-p.sigmaT * p.get_x() / p.get_mu());
}

/*----------------------------------------------------------------------------------*/
/* CAS ABSORBANT HOMOGENE SOURCE UNIFORME */
//tirage de départ de neutrons dans un matériau homogène, source uniforme
//arrivés en (x,mu).
vector_points no_scattering_homog_unif_MC(int N, point p)
{
    vector_points selection(N);
    for (int i = 0; i < N; i++)
    {
        selection.points[i] = deplacement_x(point(new_x(), p.get_mu()));
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source isotrope uniforme et matériau homogène.
double density_no_scattering_homog_unif_MC(point p, int N)
{
    vector_points selection(N);
    selection = no_scattering_homog_unif_MC(N, p);
    double dx = 1. / 100;
    double rsl = 0.;
    for (int i = 0; i < N; i++)
    {
        if (selection.points[i].get_x() > p.get_x() && selection.points[i].get_x() < p.get_x() + dx)
        {
            rsl++;
        }
    }
    //print_vector(selection);
    return rsl / (N * p.sigmaT * dx);
}

// SOURCE UNIFORME. Tirage de la densité neutronique pour mu fixé sur un maillage de nb_points de l'intervalle [0,1]. Sauve phi et la fréquence dans des fichiers .txt
vector<double> density_segment_no_scattering_homog_unif_MC(int N, int nb_points, point p)
{

    double dx = 1. / nb_points;                    //largeur des intervalles
    vector_points selection(N);                    //stocke les N tirages aléatoires
    vector<double> freq(nb_points, 0.);            //stocke la fréqence pour chaque intervalle
    selection = no_scattering_homog_unif_MC(N, p); // libre parcours depuis p
    int indice = 0;                                //initialisation

    for (int i = 0; i < selection.points.size(); i++)
    {
        indice = floor(selection.points[i].get_x() / dx); //indice de l'intervalle où est le neutron
        if (indice < nb_points)
        {
            freq[indice] += (1. / (point::sigmaT * dx)) * (1. / (float)N);
        }
    }

    // enregistrement des valeurs de phi
    ofstream fichier("Data/phi_q5_" + to_string(N) + ".txt", ios::out | ios::trunc);
    for (int i = 0; i < freq.size(); i++)
    {
        fichier << freq[i] << ",";
    }
    fichier.close();

    // enregistrement des valeurs de fréquence pour les intervalles de confiance
    ofstream fichier2("Data/freq_q5_" + to_string(N) + ".txt", ios::out | ios::trunc);
    for (int i = 0; i < freq.size(); i++)
    {
        fichier2 << freq[i] * N * point::sigmaT * dx << ",";
    }
    fichier2.close();

    return freq;
}

//calcul de la solution analytique sans scattering, cas homogène
//et source uniforme : disjonction cas sur signe de mu
double density_no_scattering_homog_unif(point p)
{
    if (p.get_mu() == 0)
    {
        return 1;
    }
    if (p.get_mu() < 0)
    {
        return (1. / p.sigmaT) * (1 - exp(-p.sigmaT * (p.get_x() - 1) / p.get_mu()));
    }
    return (1. / p.sigmaT) * (1 - exp(-p.sigmaT * p.get_x() / p.get_mu()));
}

/*----------------------------------------------------------------------------------*/
/* CAS DIFFUSANT HOMOGENE SOURCE UNIFORME */

//Calcul de la densité neutronique moyennée en mu sur un maillage de nb_points points
vector<double> density_tilda_segment_scattering_homg_unif_MC(int N, int nb_intervalles, int max_iter, double epsilon)
{

    vector_points selection(N);
    vector_points new_selection(N);
    double dx = 1. / nb_intervalles;
    vector<double> phi(nb_intervalles); // valeurs de l'estimateur
    vector<double> freq(nb_intervalles, 0.);
    vector<double> phi_old(nb_intervalles, 1.);

    //initialisations des variables
    int step = 0;
    int indice = 0;
    int stop = 0; // VA d'arret pour chaque étape et chaque tirage
    double new_x_ = 0.;
    int count = 0;

    //initalisation des positions et angles
    for (int i = 0; i < N; i++)
    {
        selection.points[i] = point(new_x(), new_mu());
    }

    while (compare_vect(phi, phi_old) > epsilon && step < max_iter && selection.points.size() > 0)
    {
        //mise à jour pour la nouvelle itération
        new_selection.points.resize(0);
        count = 0;
        phi_old = phi;

        // génération du flux n+1 fois collisionné
        for (int i = 0; i < selection.points.size(); i++)
        {
            stop = do_I_stop(point::sigmaS, point::sigmaA);
            new_x_ = deplacement_x(selection.points[i]).get_x();
            if (new_x_ >= 0 && new_x_ <= 1 && stop == 1)
            {
                new_selection.points.push_back(point(new_x_, new_mu()));
                count++;
            }
        }
        freq.resize(new_selection.points.size());
        selection = new_selection;

        //calcul de la valeur de ce flux sur chaucn des intervalles
        for (int i = 0; i < selection.points.size(); i++)
        {
            indice = floor(selection.points[i].get_x() / dx); //indice de l'intervalle où se trouve la particule
            freq[indice] += (1. / (2. * point::sigmaT * dx)) * (1. / (float)N);
        }

        //mise à jour de phi
        for (int i = 0; i < phi.size(); i++)
        {
            phi[i] += freq[i];
            freq[i] = 0;
        }
        step++;
    }
    cout << "Convergé en " << step << " itérations" << endl;

    //sauver la valeur de phi dans un fichier txt
    ofstream fichier("Data/phi_q8_" + to_string(N) + ".txt", ios::out | ios::trunc);
    for (int i = 0; i < phi.size(); i++)
    {
        fichier << phi[i] << ",";
    }
    fichier.close();
    return phi;
}

vector<double> density_tilda_segment_scattering_woodcock_unif_MC(int N, int nb_intervalles, int max_iter, double epsilon)
{

    vector_points selection(N);
    vector_points new_selection(N);
    double dx = 1. / nb_intervalles;
    vector<double> phi(nb_intervalles); // valeurs de l'estimateur
    vector<double> freq(nb_intervalles, 0.);
    vector<double> phi_old(nb_intervalles, 1.);

    //initialisations des variables
    int step = 0;
    int indice = 0;
    int stop = 0; // VA d'arret pour chaque étape et chaque tirage
    double new_x_ = 0.;
    int count = 0;

    //initalisation des positions et angles
    for (int i = 0; i < N; i++)
    {
        selection.points[i] = point(new_x(), new_mu());
    }

    while (compare_vect(phi, phi_old) > epsilon && step < max_iter && selection.points.size() > 0)
    {
        //mise à jour pour la nouvelle itération
        new_selection.points.resize(0);
        count = 0;
        phi_old = phi;

        // génération du flux n+1 fois collisionné
        for (int i = 0; i < selection.points.size(); i++)
        {
            new_x_ = deplacement_x(selection.points[i]).get_x();
            stop = do_I_stop_woodcock(point::sigmaT, point::sigmaS, sigmaT_non_cst(point(new_x_, selection.points[i].get_mu())));
            if (new_x_ >= 0 && new_x_ <= 1 && stop != 0)
            {
                if (stop == 1)
                    new_selection.points.push_back(point(new_x_, new_mu()));
                else
                    new_selection.points.push_back(point(new_x_, selection.points[i].get_mu()));
                count++;
            }
        }
        freq.resize(new_selection.points.size());
        selection = new_selection;

        //calcul de la valeur de ce flux sur chaucn des intervalles
        for (int i = 0; i < selection.points.size(); i++)
        {
            indice = floor(selection.points[i].get_x() / dx); //indice de l'intervalle où se trouve la particule
            freq[indice] += (1. / (2. * point::sigmaT * dx)) * (1. / (float)N);
        }

        //mise à jour de phi
        for (int i = 0; i < phi.size(); i++)
        {
            phi[i] += freq[i];
            freq[i] = 0;
        }
        step++;
    }
    cout << "Convergé en " << step << " itérations" << endl;

    //sauver la valeur de phi dans un fichier txt
    ofstream fichier("Data/phi_q8_" + to_string(N) + ".txt", ios::out | ios::trunc);
    for (int i = 0; i < phi.size(); i++)
    {
        fichier << phi[i] << ",";
    }
    fichier.close();
    return phi;
}

// Calcul de la densité avec discrétisation en x et en mu : demande de tirer beaucoup de points si on veut en avoir dans l'intervalle considéré. Non utilisé.
double density_scattering_homog_unif_MC(int N, point p, int max_iter, double epsilon)
{
    vector_points selection(N);
    vector_points new_selection(N);
    double dx = 1. / 100;
    double dmu = sqrt(2. / N); // discrétisation des angles
    double phi = 0.;           // valeur de l'estimateur
    double freq = 0.;
    double phi_old = 1.; //ancien phi
    int step = 0;
    int stop = 0; // VA d'arret pour chaque étape et chaque tirage
    double new_x_ = 0.;
    int count = 0;

    //initalisation
    for (int i = 0; i < N; i++)
    {
        selection.points[i] = point(new_x(), new_mu());
    }

    while (abs(phi - phi_old) > epsilon && step < max_iter && selection.points.size() > 0)
    {
        new_selection.points.resize(0);
        count = 0;
        for (int i = 0; i < selection.points.size(); i++)
        {
            stop = do_I_stop(point::sigmaS, point::sigmaA);
            new_x_ = deplacement_x(selection.points[i]).get_x();

            if (stop == 1 && new_x_ >= 0 && new_x_ <= 1)
            {
                new_selection.points.push_back(point(new_x_, new_mu()));
                count++;
            }
        }
        selection = new_selection;

        for (int i = 0; i < selection.points.size(); i++)
        {
            if (selection.points[i].get_x() > p.get_x() - dx / 2 && selection.points[i].get_x() < p.get_x() + dx / 2 && selection.points[i].get_mu() > p.get_mu() && selection.points[i].get_mu() < p.get_mu() + dmu)
            {
                freq++;
            }
        }

        phi_old = phi;
        phi += (1. / (2. * point::sigmaT * dx * dmu)) * (freq / (float)N);
        freq = 0.;
        step++;
    }
    cout << "Convergé en " << step << " itérations" << endl;
    return phi;
}
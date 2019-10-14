#include "aleat_tools.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <typeinfo>

using namespace std;

//tirage d'un x uniforme entre 0 et 1
double new_x()
{
    return (static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
}

//tirage d'une variable aléatoire pour le parcours entre x_n-1 et x_n
double parcours_x(point p)
{
    double y = new_x(); // loi uniforme sur (0,1)
    return -(p.get_mu() / p.sigmaT) * log(1 - y);
}

//nouveau point p après déplacement
point deplacement_x(point p)
{
    double new_x = p.get_x() + parcours_x(p);
    point p2(new_x, p.get_mu());
    return p2;
}

//tirage d'une variable aléatoire pour la continuation de la marche 0 arret, 1 continuer
double do_I_stop(double sigmaS, double sigmaA)
{
    double y = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (sigmaA + sigmaS);
    if (y < sigmaA)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

//tirage d'un cos(angle) uniforme entre -1 et 1
double new_mu()
{
    return 2 * (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) - 1.;
}

//maximum d'un vecteur de double
double vector_max(vector<double> v){
    double maxi = 0.;
    for(int i=0;i<v.size();i++){
        if(maxi<=v[i]) {
            maxi=v[i];
        }
    }
    return maxi;
}

//norme infini de la différence de deux vecteurs
double compare_vect(vector<double> v1, vector<double> v2){
    if (v1.size()!= v2.size()) {return 0;}
    else {
        double maxi=0.;
        for(int i=0;i<v1.size();i++){
            if(maxi<=abs(v1[i]-v2[i])) {
            maxi=abs(v1[i]-v2[i]);
            }
        }
        return maxi;
    }
}
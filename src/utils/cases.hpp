#ifndef CASES_HPP_INCLUDED
#define CASES_HPP_INCLUDED

#include <vector>
#include <math.h>
#include <iostream>
#include "points.hpp"
#include "aleat_tools.hpp"

class source_point
{
public:
    point change_x(point p);
};

point source_point::change_x(point p)
{
    return parcours_x_p(point(0, p.get_mu()));
}

class source_uniforme
{
public:
    point change_x(point p);
};

point source_uniforme::change_x(point p)
{
    return deplacement_x(point(new_x(), p.get_mu()));
}

template <typename SOURCE>
double density_no_scattering_MC(double mu, double x, int N)
{
    double dx = 1. / 100;
    double rsl = 0.;
    double pt = 0.;
    SOURCE source;
    for (int i = 0; i < N; i++)
    {
        pt = source.change_x(point(x, mu)).get_x();
        if (pt > x - dx / 2 && pt < dx / 2 + x)
            rsl++;
    }
    return rsl / (N * point::sigmaT * dx);
}

class absorption
{
public:
    int stop_test(double sigmaS, double sigmaA, int iter) { return (iter==0)?1:0; }
    template <typename SOURCE>
    point new_point_(point p) {return p;}
};

class scattering
{
public:
    int stop_test(double sigmaS, double sigmaA, int iter) { return do_I_stop(sigmaS, sigmaA); }
    template <typename SOURCE>
    point new_point_(point p) {
        SOURCE source ;
        return source.change_x(p);}
};

template <typename SOURCE, typename PHYSIQUE>
vector<double> density_MC(int N, int nb_intervalles, int max_iter)
{
    double dx = 1. / nb_intervalles;
    vector<double> freq(nb_intervalles, 0.);
    SOURCE source;
    PHYSIQUE physique;

    //initialisations des variables
    int step = 0;
    int indice = 0;
    point new_point;
    int stop = 1;

    for (int i = 0; i < N; i++)
    {
        new_point = point(new_x(), new_mu());

        while (step < max_iter && stop != 0)
        {
            new_point =  physique.new_point_(new_point) ;

            if (new_point.get_x() >= 0 && new_point.get_x() <= 1 && physique.stop_test(point::sigmaS, point::sigmaA, step) == 1)
            {
                indice = floor(new_point.get_x() / dx);
                freq[indice] += (1. / (point::sigmaT * dx)) * (1. / (float)N);
            }
            else
                stop = 0;
            step++;
        }
        step = 0;
        stop = 1;
    }
    return freq;
}

#endif // CASES_HPP_INCLUDED
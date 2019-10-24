#ifndef DETERMINISTE_INCLUDED
#define DETERMINISTE_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include "points.hpp"

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

vector<double> IS_iteration(int Nx, double mu, double source, vector<double> Q,vector<double> sigmaT);
vector<double> IS(int Nx, int Nmu, double epsilon, int iter_max, vector<double> S, vector<double> sigmaT);

#endif // DETERMINISTE_INCLUDED
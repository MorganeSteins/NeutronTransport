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

vector<double> solveur(int Nx, double mu, double source);


#endif // DETERMINISTE_INCLUDED
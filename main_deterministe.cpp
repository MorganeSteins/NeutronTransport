#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include "points.hpp"
#include "aleat_tools.hpp"
#include "MC.hpp"
#include "deterministe.hpp"

int main() {

    double mu=1.;
    int Nx=50;
    vector<double> phi = solveur(Nx,mu,1./mu);
    cout<<phi[Nx]<<endl;
    return 0;
}
#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "src/utils/points.hpp"
#include "src/utils/aleat_tools.hpp"
#include "src/mc/MC.hpp"
#include "src/mc/Mc_iter.hpp"
// #include "matplotlibcpp.h"
#include "src/deterministe/deterministe.hpp"
// #include "cases.hpp"

#include <chrono>
#include <unistd.h>

using namespace std;
// namespace plt = matplotlibcpp;

int main(int argc, char *argv[])
{
    srand(static_cast<unsigned>(time(0)));

    // reading inputs
    double x = 1, mu = 1;
    int N = 1000000;
    int max_iter = 50;
    string type_ = "density";
    if (argc == 1)
    {
        printf("You did not mention the type of requested value, going with density");
    }
    else
    {
        type_ = argv[1];
    }
    if (argc < 5)
    {
        printf("Not enough inputs, going with default values N=10^6, max_iter=50\n");
    }
    else
    {
        N = atof(argv[2]);
        max_iter = atof(argv[3]);
        if (type_.compare("point") == 0 && argc > 3)
        {
            x = atof(argv[4]);
            mu = atof(argv[5]);
        }
        else
        {
            printf("You requested point value evaluation but did not mention the point, giving answer for x=1 and mu=1\n");
        }
    }

    // // QUESTION 4
    // // Compute density at point (x,mu)
    cout << "Q4 : Compute phi at (" << x << "," << mu << ") with N=" << N << endl;
    cout << "densité par MC " << density_no_scattering_homog_point_MC(point(x, mu), N) << endl;
    cout << "densité réelle " << density_no_scattering_homog_point(point(x, mu)) << endl;

    point p(x, mu);
    auto start = chrono::steady_clock::now();
    auto data1 = density_no_scattering_homog_point_MC(p, N);
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    start = chrono::steady_clock::now();
    auto data2 = density_no_scattering_homog_point_MC_iter(p, N);
    end = chrono::steady_clock::now();
    cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";
    cout << "Error comparison " << data2 - data1 << "\n\n";

    cout << "\n ----TEST 2---- \n";
    start = chrono::steady_clock::now();
    data1 = density_tilda_no_scattering_homog_point_MC(p, N);
    end = chrono::steady_clock::now();
    cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    start = chrono::steady_clock::now();
    data2 = density_tilda_no_scattering_homog_point_MC_iter(p, N);
    end = chrono::steady_clock::now();
    cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n";

    start = chrono::steady_clock::now();
    auto data3 = density_tilda_no_scattering_homog_point_MC_iter(x, N);
    end = chrono::steady_clock::now();
    cout << "Elapsed time " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms \n\n";
    cout << "Error comparison " << data2 - data1 << "\n";
    cout << "Error comparison " << data3 - data1 << "\n\n";

    // // QUESTION 5
    // // Compute density at (x,mu)
    cout << "Q5 density at point (" << x << "," << mu << ") with N=" << N << endl;
    cout << "MC value " << density_no_scattering_homog_unif_MC(point(x, mu), N) << endl;
    cout << "Ground truth" << density_no_scattering_homog_unif(point(x, mu)) << endl;

    // //calcul de la densité en nb_points intervalles et sauvée dans un fichier txt pour les plots
    int nb_points = 25;
    cout << "N = " << N << " and  there is " << nb_points << " intervalles." << endl
         << endl;
    density_segment_no_scattering_homog_unif_MC(N, nb_points, p);

    // QUESTION 8
    //Only averaged over mu values of phis

    double value_point = density_scattering_homog_unif_MC(1e7, p, max_iter);
    vector<double> value_h = density_tilda_segment_scattering_homg_unif_MC(N, nb_points, max_iter);
    vector<double> value_wd = density_tilda_segment_scattering_woodcock_unif_MC(N, nb_points, max_iter);
    cout << "Q8 value at 0 is " << value_h[0] << " and value with woodcock at 0 is " << value_wd[0] << endl;
    cout << "Complete values are stored in files for python plot precessing" << endl
         << endl;

    return 0;
}

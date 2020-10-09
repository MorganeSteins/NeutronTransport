#ifndef POINTS_HPP_INCLUDED
#define POINTS_HPP_INCLUDED

#include <vector>
#include <iostream>

using namespace std;

class point
{
    double x;
    double mu;

public:
    point()
    {
        x = 0;
        mu = 0;
    }
    point(double newx, double newmu)
    {
        x = newx;
        mu = newmu;
    }
    double get_x() const { return x; }
    double get_mu() const { return mu; }
    void set_x(double x1) { x = x1; }
    void set_mu(double mu1) { mu = mu1; }

    static double sigmaS;
    static double sigmaA;
    static double sigmaT;
};

class vector_points
{
public:
    vector<point> points;

    vector_points() { points = vector<point>(0); };
    vector_points(int dimension);
    vector_points(vector_points &vector_to_copy);

    vector<point> get_points() const { return points; }
    void add_point(point p);
    void iteration_MC();
};

void print_vector(vector_points &v);

double sigmaT_non_cst(point);

#endif // POINTS_HPP_INCLUDED
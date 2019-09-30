#ifndef POINTS_HPP_INCLUDED
#define POINTS_HPP_INCLUDED

#include <vector>
#include <iostream>

using namespace std;

class point {
    double x;
    double mu;

    public :
    point() {x=0;mu=0;}
    point(double newx, double newmu) {x=newx;mu=newmu;}
    double get_x() const {return x;}
    double get_mu() const {return mu;}
    double set_x(double x1) {x=x1;}
    double set_mu(double mu1) {mu=mu1;}

    static double sigmaS;
    static double sigmaA;
    static double sigmaT;
};

class vector_points {
public :
    int dim;
    vector<point> points ;

    vector_points() {dim=0;points= vector<point>(0) ;};
    vector_points(int dimension) ;
    vector_points(vector_points &vector_to_copy) ;

    int dimension() const {return dim;}
    vector<point> get_points() const {return points;}
    void set_points(int i, point p)  { if (i>dim-1){dim=i+1;points.resize(dim, point());} points[i]=p;}
    void add_point(point p);
    void iteration_MC();

};

void print_vector(vector_points &v);

#endif // POINTS_HPP_INCLUDED
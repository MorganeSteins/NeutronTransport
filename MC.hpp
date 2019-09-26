#ifndef MC_HPP_INCLUDED
#define MC_HPP_INCLUDED

#include <vector>

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
};

class vector_points {
    int dim;
    vector<point> points ;

public :
    vector_points() {dim=0;points= vector<point>(0) ;};
    vector_points(int dimension) {dim=dimension;points= vector<point>(dimension) ;};
    vector_points(vector_points &list_to_copy) ;
    //~vector_points() {delete points;}

    int dimension() const {return dim;}
    vector<point> get_points() const {return points;}
    vector<point> set_points(int i, point p)  { if (i<dim-1){points[i]=p;}}
};

//faire une classe liste_points pour stocker toutes les contributions

double parcourt_x(point p, double sigmat);

#endif // MC_HPP_INCLUDED
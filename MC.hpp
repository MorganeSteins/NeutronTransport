#ifndef MC_HPP_INCLUDED
#define MC_HPP_INCLUDED

class point {
    double x;
    double mu;

    public :
    point() {x=0;mu=0;}
    point(double newx, double newmu) {x=newx;mu=newmu;}
    double get_x() {return x;}
    double get_mu() {return mu;}

};

//faire une classe liste_points pour stocker toutes les contributions

double parcourt_x(point p, double sigmat);

#endif // MC_HPP_INCLUDED
#include "MC.hpp"
#include <math.h>
#include<iostream>
#include <vector>
#include <cstdlib>

using namespace std;

//initialisation des variabales statiques
double vector_points::sigmaA = 0.1;
double vector_points::sigmaS = 0.9;

//créateur par copie
vector_points::vector_points(vector_points &vector_to_copy) {
    dim = vector_to_copy.dimension();
    vector<point> points(dim);
    for (int i=0;i<dim;i++) {this->points[i]=vector_to_copy.points[i];}

}

//createur d'un vecteur de points initialisés à (0,0)
vector_points::vector_points(int dimension) {
    //if (dimension<0) {throw length_error{"Vector::Vector"}}; //jsp comment faire
    dim = dimension;
    points = vector<point>(dimension) ;
    for (int i=0;i<dimension;i++){
        points[i] = point(0.,0.);
    }
}

//Push_Back nouveau point à la fin du vecteur et change la dimension
void vector_points::add_point(point p){
    dim++;
    points.push_back(p);
}

//affiche le vecteur de points 
void print_vector(vector_points &v){
    int dim = v.dimension();
    cout<<3<<endl;
    vector<point> points=v.get_points();
    cout<<"[";
    for (int i=0;i<dim-1;i++){
        cout<<"("<<points[i].get_x()<<" , "<<points[i].get_mu()<<");"<<endl;
    }
    cout<<"("<<points[dim-1].get_x()<<" , "<<points[dim-1].get_mu()<<")";
    cout<<"]"<<endl;
}

//FOCNTIONS AUXILIAIRES

//tirage d'une variable aléatoire pour le parcourt entre x_n-1 et x_n
double parcourt_x(point p, double sigmaT) {
    double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); // loi uniforme sur (0,1)
    return -(p.get_mu()/sigmaT)*log(1-y);
}

//nouveau point p après déplacement
point deplacement_x(point p, double sigmaT){
    double new_x =  p.get_x()+p.get_mu()*parcourt_x(p, sigmaT);
    point p2;
    p2.set_mu(p.get_mu());
    p2.set_x(new_x);
    return p2;
}

//tirage d'une variable aléatoire pour la continuation de la marche 0 arret, 1 continuer
double do_I_stop(double sigmaS, double sigmaA) {
    double y = (static_cast <double> (rand()) / static_cast <double> (RAND_MAX))*(sigmaA+sigmaS); 
    if (y<sigmaA) {return 0;}
    else {return 1;}
}

//tirage d'un cos(angle) uniforme entre -1 et 1
double new_mu() {
    return 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX))-1.;
}




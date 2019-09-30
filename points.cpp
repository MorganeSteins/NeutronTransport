#include "points.hpp"
#include <math.h>
#include<iostream>
#include <vector>
#include <cstdlib>

using namespace std;

//initialisation des variabales statiques
double point::sigmaA = 1.;
double point::sigmaS = 0.;
double point::sigmaT = sigmaA+sigmaS;

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
    vector<point> points=v.get_points();
    cout<<"[";
    for (int i=0;i<dim-1;i++){
        cout<<"("<<points[i].get_x()<<" , "<<points[i].get_mu()<<");"<<endl;
    }
    cout<<"("<<points[dim-1].get_x()<<" , "<<points[dim-1].get_mu()<<")";
    cout<<"]"<<endl;
}


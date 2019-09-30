#include <iostream>
#include <math.h>
#include "points.hpp"
#include "aleat_tools.hpp"
#include <vector>
#include <fstream>
#include <time.h>

using namespace std;

int main() {
    srand (static_cast <unsigned> (time(0)));
    double sigmat=1;
    point a(1,1);
    cout<<a.get_x()<<" , "<<a.get_mu()<<endl;


    vector_points b(1);
    cout<<b.dimension()<<endl;
    cout<<"mu vaut "<<b.get_points()[0].get_mu()<<endl;

    b.set_points(0,a);
    cout<<"mu vaut "<<b.get_points()[0].get_mu()<<" et le x "<<b.get_points()[0].get_x()<<endl;;
    
    b.add_point(point(2,2));
    cout<<"b est de dimension "<<b.dimension()<<endl;
    cout<<"le second mu de b vaut "<<b.get_points()[1].get_mu()<<" et le x "<<b.get_points()[1].get_x()<<endl;
    
    cout<<"sigmaT vaut"<<point::sigmaT<<endl;
 
    return 0;
}
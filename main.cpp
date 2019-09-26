#include <iostream>
#include <math.h>
#include "MC.hpp"
#include <vector>
#include <fstream>
#include <time.h>

using namespace std;

int main() {
    double sigmat;
    point a(1,1);
    cout<<a.get_x()<<" , "<<a.get_mu()<<endl;
    vector_points b(1);
    b.set_points(0,a);
    cout<<b.get_points()[0].get_mu()<<endl;
    return 0;
}
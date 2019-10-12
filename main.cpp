#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include "points.hpp"
#include "aleat_tools.hpp"
#include "MC.hpp"

using namespace std;

int main() {
    srand (static_cast <unsigned> (time(0)));
    
    //print_vector(sol);
    double x = 0.5;
    double mu = 0.1;
    double integrale = 0.;
    cout<<"On travaille à x="<<x<<" et mu="<<mu<<" avec sigmaT "<<point(0,0).sigmaT<<endl;

    // test question 4
    /*double density = density_no_scattering_homog_point_MC(point(x,mu), 1000000);
    cout<<"densité par MC "<<density<<endl;
    cout<<"densité réelle "<<density_no_scattering_homog_point(point(x,mu))<<endl<<endl;
    */
   //save_error(1000000,point(x,mu));
   //save_points(1000000,25,mu, "points_q4.txt");
   string filename_;
   int nb_points = 25;
   int N = 1000000;
   for (int i=1;i<=N;i=i*10){
        filename_ = "Data/points_mu01_q5_"+to_string(i)+".txt";
        cout<<"Je calcule pour "<<i<<" points"<<endl;
        save_points(i,nb_points,mu,filename_);
   }


    // test question 5
    /*double density_2 = density_no_scattering_homog_unif_MC(point(x,mu), 1000000);
    cout<<"densité par MC "<<density_2<<endl;
    cout<<"densité réelle "<<density_no_scattering_homog_unif(point(x,mu))<<endl<<endl;
    */

    // tests vecteurs
    /*vector_points v1(2);
    v1.points[0] = point(1,1);
    v1.points[1] = point(1,2);
    vector_points v2(1);
    v2.points[0] = point(1,0);
    cout<<v1.points.size()<<endl;
    v1=v2;
    v1.points.resize(v2.points.size());
    cout<<v1.points.size()<<endl;
    cout<<v1.points[1].get_mu()<<endl;

    vector<int> v(1000);
    generate(v.begin(), v.end(), std::rand);
    cout<<v[0]<<endl;*/

    // test question 8
    //double density_scattering = density_tilda_scattering_homg_unif_MC(1000000,x,10000, 0.0001);
    //cout<<" Densité MC "<<density_scattering<<endl;
    /*vector_points vp(2);
    vp.points[0] = point(0.5,0.5);
    vp.points[1] = point(0.2,0.2);
    vector_points sol(0);*/
    

    //for (int i=0;i<1000;i++) {integrale += (15./1000.)*density_no_scattering_homog_unif(point((float) i*15/1000,mu));}
    //cout<<integrale<<endl;
    //for (int i=0;i<11;i++) {cout<<(float)i/10<<" vaut "<<density_no_scattering_homog_unif_MC(point(float(i)/10,mu),10000)<<endl;}
    /* test du flux non collisionné : intégrale somme à 1*/

    // somme des valeurs pour approcher l'intégrale MC : doit faire 1. De -7 à 7.
    /*for (int i=0;i<1000;i++) {integrale += (1./1000.)*density_no_scattering_homog_unif(point((float) i*1./1000,mu));}
    cout<<integrale<<endl;*/
    return 0;
}

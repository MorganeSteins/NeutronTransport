#include "MC.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

/* CAS ABSORBANT HOMOGENE SOURCE PONCTUELLE */
//tirage de départ de neutrons dans un matériau homogène, source ponctuelle en 0.
vector_points no_scattering_homog_point_MC(int N, double mu){
    vector_points selection(N);
    for (int i=0; i<N; i++){
        selection.points[i] = point(parcours_x(point(0,mu)),mu);
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source ponctuelle isotrope en 0 et matériau homogène.
double density_no_scattering_homog_point_MC(point p, int N){
    double dx=sqrt(1./N);
    vector_points selection(N) ;
    selection = no_scattering_homog_point_MC(N,p.get_mu());
    double rsl = 0.;
    for (int i=0;i<N;i++){
        if (selection.points[i].get_x()>p.get_x() && selection.points[i].get_x()<dx+p.get_x()) {
            rsl++;
        }
    }
    return rsl/(N*p.sigmaT*dx);
}

// calcul de la densité neutronique au point x moyennée en mu avec N tirages
// avec une source ponctuelle isotrope en 0 et matériau homogène.
double density_tilda_no_scattering_homog_point_MC(point p, int N){
    double dx=sqrt(1./N);
    vector_points selection(N) ;
    selection = no_scattering_homog_point_MC(N,new_mu());
    double rsl = 0.;
    for (int i=0;i<N;i++){
        if (selection.points[i].get_x()>p.get_x() && selection.points[i].get_x()<dx+p.get_x()) {
            rsl++;
        }
    }
    return rsl/(2*N*p.sigmaT*dx);
}

//calcul de la solution analytique sans scattering, cas homogène 
//et source ponctuelle
double density_no_scattering_homog_point(point p){
    if (p.get_mu()==0) {if (p.get_x()==0){return 1;} return 0;}
    return (1/p.get_mu())*exp(-p.sigmaT*p.get_x()/p.get_mu());}



/* CAS ABSORBANT HOMOGENE SOURCE UNIFORME */ 
//tirage de départ de neutrons dans un matériau homogène, source uniforme 
//arrivés en (x,mu).
vector_points no_scattering_homog_unif_MC(int N, point p){
    vector_points selection(N);
    for (int i=0; i<N; i++){
        selection.points[i] = deplacement_x(point(new_x(), p.get_mu()));
    }
    return selection;
}

// calcul de la densité neutronique au point x direction mu avec N tirages
// avec une source isotrope uniforme et matériau homogène.
double density_no_scattering_homog_unif_MC(point p, int N){
    vector_points selection(N) ;
    selection = no_scattering_homog_unif_MC(N,p);
    double dx = sqrt(1./N);
    double rsl = 0.;
    for (int i=0;i<N;i++){
        if (selection.points[i].get_x()>p.get_x() && selection.points[i].get_x()<p.get_x()+dx) {
            rsl++;
        }
    }
    //print_vector(selection);
    return rsl/(N*p.sigmaT*dx);
}

//calcul de la solution analytique sans scattering, cas homogène 
//et source uniforme : disjonction cas sur signe de mu 
double density_no_scattering_homog_unif(point p){
    if (p.get_mu()==0) {return 1;}
    if (p.get_mu()<0) {return (1./ p.sigmaT)*(1-exp(-p.sigmaT*(p.get_x()-1)/p.get_mu()));}
    return (1./ p.sigmaT)*(1-exp(-p.sigmaT*p.get_x()/p.get_mu()));}



/* CAS DIFFUSANT HOMOGENE SOURCE UNIFORME */
/* Génération des nouveaux points (x_n+1,mu_n) à partir de (x_n,mu_n) 
 Ne conserve que ceux qui restent dans l'intervalle [0,1] */
/*void move_scattering(vector_points start,vector_points final){
    cout<<"je commence "<<endl;
    double new_x_;
    int stop=0;
    int count=0;
    cout<<"La taille du vecteur max est "<<start.points.size()<<endl;
    for (int i=0;i<start.points.size();i++){
        new_x_ = deplacement_x(start.points[i]).get_x();
        stop = do_I_stop(point::sigmaS,point::sigmaA);
        if (new_x_>=0 && new_x_<=1 && stop==1){
            final.points.push_back(point(new_x_,new_mu()));
            count ++;
        }
    }
    cout<< " I y a eu "<<count<<" points non nuls"<<endl;
    final.points.resize(count);
}*/


double density_tilda_scattering_homg_unif_MC(int N, double x, int max_iter, double epsilon ){
    vector_points selection(N);
    vector_points new_selection(N);
    double dx = sqrt(1./N);
    double phi=0.; // valeur de l'estimateur
    double freq = 0.;
    double phi_old = 1.; //ancien phi
    int step = 0;
    int stop=0; // VA d'arret pour chaque étape et chaque tirage
    double new_x_=0.;
    int count=0;

    //initalisation 
    for (int i=0;i<N;i++) {
        selection.points[i] = point(new_x(),new_mu());
    }

    while (abs(phi-phi_old)>epsilon && step<max_iter && selection.points.size()>0){
        //cout<<"Phi a l'itération "<<step<<" vaut "<<phi<<endl;
        new_selection.points.resize(0);
        count = 0;
        for (int i=0;i<selection.points.size();i++){
            stop = do_I_stop(point::sigmaS,point::sigmaA);
            new_x_ = deplacement_x(selection.points[i]).get_x();
            //cout<<" Le nouveau x vaut "<<new_x_<<endl;
            if (new_x_>=0 && new_x_<=1 && stop==1){
                new_selection.points.push_back(point(new_x_,new_mu()));
                count ++;
            }
        }
        selection = new_selection;
        for (int i=0;i<selection.points.size();i++) {
            if (selection.points[i].get_x()>x && selection.points[i].get_x()<x+dx){
                freq ++;
            }
        }
        phi_old = phi;
        phi += (1./(2.*point::sigmaT*dx))*(freq/(float)N);
        //cout<<phi_old<<" et "<<phi<<endl<<endl;
        freq = 0.;
        step++;
    }
    cout<<"Convergé en "<<step<<" itérations"<<endl;
    return phi;
}

double density_scattering_homg_unif_MC(int N, point p, int max_iter, double epsilon ){
    vector_points selection(N);
    vector_points new_selection(N);
    double dx = 1./100;
    double dmu = sqrt(2./N); // discrétisation des angles
    double phi=0.; // valeur de l'estimateur
    double freq = 0.;
    double phi_old = 1.; //ancien phi
    int step = 0;
    int stop=0; // VA d'arret pour chaque étape et chaque tirage
    double new_x_=0.;
    int count=0;

    //initalisation 
    for (int i=0;i<N;i++) {
        selection.points[i] = point(new_x(),new_mu());
    }

    while (abs(phi-phi_old)>epsilon && step<max_iter && selection.points.size()>0){
        //cout<<"Phi a l'itération "<<step<<" vaut "<<phi<<endl;
        new_selection.points.resize(0);
        count = 0;
        for (int i=0;i<selection.points.size();i++){
            stop = do_I_stop(point::sigmaS,point::sigmaA);
            new_x_ = deplacement_x(selection.points[i]).get_x();
            //cout<<" Le nouveau x vaut "<<new_x_<<endl;
            if (stop==1 && new_x_>=0 && new_x_<=1){
                new_selection.points.push_back(point(new_x_,new_mu()));
                count ++;
            }
        }
        selection = new_selection;
        for (int i=0;i<selection.points.size();i++) {
            if (selection.points[i].get_x()>p.get_x() && selection.points[i].get_x()<p.get_x()+dx && selection.points[i].get_mu()>p.get_mu() && selection.points[i].get_mu()<p.get_mu()+dmu){
                freq ++;
            }
        }
        cout<<"Iteration "<<step<<" il reste "<<selection.points.size()<<" il ya eu "<<freq<<endl;
        phi_old = phi;
        phi += (1./(2.*point::sigmaT*dx*dmu))*(freq/(float)N);
        freq = 0.;
        step++;
        
    }
    cout<<"Convergé en "<<step<<" itérations"<<endl;
    return phi;
}


void save_error(int N_max,point p)
{
    ofstream fichier("erreurs_q4.txt", ios::out | ios::trunc);
    double phi_exact = density_no_scattering_homog_point(p);
    double err=0.;
    if (fichier) // si l'ouverture a reussi
    {
        //on ecrit dans le fichier
        cout << "Calcul des erreurs pour différents N jusqu'à "<<endl<<N_max << endl;
        //on écrit les erreurs
        for (int i = 1; i <= N_max; i=i*10){
            cout<<"i="<<i<<" et "<<2*i<<endl;
            err=0.;
            for (int j=0;j<5;j++){err+=abs(density_no_scattering_homog_point_MC(p,i)-phi_exact);}
            fichier << err/5<<",";
            err=0.;
            for (int j=0;j<5;j++){err+=abs(density_no_scattering_homog_point_MC(p,2*i)-phi_exact);}
            fichier << err/5<<",";
        }
        fichier.close(); // on referme le fichier
    }
    else {cerr << "Erreur a l'ouverture !" << endl;}
}

void save_points(int N,int nb_points,double mu, string filename_)
{
    double dx = 1./nb_points;
    double val=0.;
    ofstream fichier(filename_, ios::out | ios::trunc);
    if (fichier) // si l'ouverture a reussi
    {
        //on écrit les erreurs
        for (int i =0; i <nb_points; i++){
            val=0.;
            for (int j=0;j<5;j++){val+=density_no_scattering_homog_unif_MC(point(i*dx,mu),N);}
            fichier << val/5<<",";
        }
        fichier.close(); // on referme le fichier
    }
    else {cerr << "Erreur a l'ouverture !" << endl;}
}
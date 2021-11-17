#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

#define nx 150
#define ny 100
double delta=0.1;
double V1=10;
double V2=0;
double eps=1;
double TOL=pow(10,-8);
double xmax=delta*nx;
double ymax=delta*ny;
double sigx=0.1*15;
double sigy=0.1*10;
double density[nx+1][ny+1];


void egdge_conditions(double V[][ny+1]){
    for(int i=0;i<=nx;i++){
        V[i][0]=V1;
        V[i][ny]= V2;
    }
}

void potential_density(){

    for(int i=0;i<=nx;i++){
        for(int j= 0;j<=ny;j++){
		    	density[i][j]=exp(-pow(delta*i-0.35*xmax,2)/(pow(sigx,2)) - pow(delta*j-0.5*ymax,2)/(pow(sigy,2)));
            density[i][j]+=(-1.0)*exp(-pow(delta*i-0.65*xmax,2)/(pow(sigx,2))-pow(delta*j-0.5*ymax,2)/(pow(sigy,2)));
            }
        }
}

double is_termination_condition(double V[][ny+1]){
    double S=0.0;
    for(int i=0;i<nx;i++) 
		for(int j=0;j<ny;j++) 
				S+=pow(delta,2)*(0.5*pow((V[i+1][j]-V[i][j])/delta,2)+0.5*pow((V[i][j+1]-V[i][j])/delta,2)-density[i][j]*V[i][j]);    	
    return S;
}


double solution_error(double V[][ny+1], int i, int j){
    return ((V[i+1][j]-2*V[i][j]+V[i-1][j])/(pow(delta,2))+(V[i][j+1]-2*V[i][j]+V[i][j-1])/(pow(delta,2))+density[i][j]/eps);
}

void global_relaxation(double omega, ofstream &integral, ofstream &map, ofstream &error){

    double S_curr,S_prev=0.;
    double Vn[nx+1][ny+1];
    double Vs[nx+1][ny+1];
    potential_density();

     for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
            Vn[i][j]=0.0; 
        }
    }
    egdge_conditions(Vn);

    for(int i=0;i<=nx;i++){
            for(int j=0;j<=ny;j++){
                Vs[i][j] =Vn[i][j];
            }
    }
    S_prev = is_termination_condition(Vn);
    int iter=0;
    do{
      iter+=1;
      for(int i=1;i<nx;i++)
        for(int j=1;j<ny;j++)
            Vn[i][j]=0.25*(Vs[i+1][j]+Vs[i-1][j]+Vs[i][j+1]+Vs[i][j-1]+pow(delta,2)/eps*density[i][j]);

      for(int j=1;j<ny;j++){
        Vn[0][j]=Vn[1][j];
        Vn[nx][j]=Vn[nx-1][j];
    }

      for(int i=0;i<=nx;i++)
        for(int j=0;j<=ny;j++)
            Vs[i][j]=(1-omega)*Vs[i][j]+omega*Vn[i][j];
      S_prev=S_curr;
      S_curr=is_termination_condition(Vs);        
      integral<<" "<<iter<<" "<<S_curr<<"\n";
    }while(fabs((S_curr-S_prev)/S_prev)>TOL);
    cout<<"GLOBAL Iteracja dla omega = "<<omega<<" "<<iter<<"\n";



    for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
           map<<delta*i<<" "<<delta*j<<" "<<Vs[i][j]<<"\n";
            if(i>0&&j>0&&i<nx&&j<ny)
               error<<delta*i<<" "<<delta*j<<" "<<solution_error(Vs, i, j)<<"\n";
             else
                error<<delta*i<<" "<<delta*j<<" "<< 0.0<<"\n";
        }
        map<<"\n";
        error<<"\n";
  }

  

}

void lokal_relaxation(double omega, ofstream &integral){

    double S_prev,S_curr=0;
    double V[nx+1][ny+1];
     for(int i=0;i <= nx; i++){
        for(int j=0;j<=ny;j++){
            V[i][j] = 0.0;
        }
    }
    egdge_conditions(V);
    S_prev = is_termination_condition(V);
    int iter = 0;
    do{
        iter+=1;
        for(int i=1;i<nx;i++)
            for(int j=1;j<ny;j++)
                V[i][j]= (1-omega)*V[i][j]+(omega/4.0)*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]+pow(delta,2)/eps*density[i][j]);

       for(int j=1;j<ny;j++){
          V[0][j]=V[1][j];
          V[nx][j]=V[nx - 1][j]; 
       }       
        S_prev=S_curr;
        S_curr=is_termination_condition(V);  
        integral<<" "<<iter<<" "<<S_curr<<"\n";
      }while(fabs((S_curr-S_prev)/S_prev)>TOL);
      cout<<"LOKAL Iteracja dla omega = "<<omega<<" "<<iter<<"\n";

}


int main(){

    double omega1[]={1.0,0.6};
    double omega2[]={1.0, 1.4, 1.8, 1.9};

    ofstream integral;
    ofstream map;
    ofstream integral1;
    ofstream map1;
    ofstream error1;
    ofstream error;
    ofstream integral2;
    ofstream integral3;
    ofstream integral4;
    ofstream integral5;

    integral.open("global_1.0.dat");
    map.open("map_1.0.dat");
    error.open("error_1.0.dat");
    integral1.open("global_0.6.dat");
    map1.open("map_0.6.dat");
    error1.open("error_0.6.dat");
    integral2.open("local_1.0.dat");
    integral3.open("local_1.4.dat");
    integral4.open("local_1.8.dat");
    integral5.open("local_1.9.dat");


    global_relaxation(omega1[0], integral, map, error);
    global_relaxation(omega1[1], integral1, map1, error1);
    
    lokal_relaxation(omega2[0], integral2);
    lokal_relaxation(omega2[1], integral3);
    lokal_relaxation(omega2[2], integral4);
    lokal_relaxation(omega2[3], integral5);
    
    integral.close();
    map.close();
    integral1.close();
    map1.close();
    error1.close();
    error.close();
    integral2.close();
    integral3.close();
    integral4.close();
    integral5.close();

    return 0;
}
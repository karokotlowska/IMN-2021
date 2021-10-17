
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

const double bet = 0.001;
const int N = 500;
const double gama = 0.1;
const int tmax = 100;
const double dt = 0.1;
const int iteracje_max = 20;
const double u0 = 1.;
const double TOL = 10e-6;
double alfa=bet*N-gama;

void pickardMethod(){
  std::ofstream f1;
  f1.open("pickard.txt");
  double u_mi=u0;
  double u_mi_1=0.;
  double un=u0; 

  for(double i=0.;i<tmax/dt;i++){
    for(int m =0;m<=iteracje_max;m++){
    u_mi_1=un+(dt/2)*((alfa*un-bet*un*un)+ (alfa*u_mi - bet*u_mi*u_mi));
      if(abs(u_mi_1 - u_mi)< TOL) break;
      else{
        u_mi = u_mi_1;
      }
    }
    un=u_mi_1;
    u_mi=u_mi_1;

    f1<<i*dt<<" "<<un<<" "<<N-un<<"\n";
  }

  f1.close();
}


void NewtonMethod(){
  std::ofstream f1;
  f1.open("newton.txt");
  double u_mi=u0;
  double u_mi_1=0.;
  double un=u0; 

  for(double i=0.;i<tmax/dt;i++){
    for(int m =0;m<=iteracje_max;m++){
    u_mi_1=u_mi-(u_mi-un-(dt/2)*((alfa*un-bet*un*un)+ (alfa*u_mi - bet*u_mi*u_mi)))/(1-(dt/2)*(alfa-2*bet*u_mi));
      if(abs(u_mi_1 - u_mi)< TOL) break;
      else{
        u_mi = u_mi_1;
      }
    }
    un=u_mi_1;
    u_mi=u_mi_1;

    f1<<i*dt<<" "<<un<<" "<<N-un<<"\n";
  }

  f1.close();
}

double m_11(double U1, double a11) {
    return (1 - dt * a11 * ((bet * N - gama) - 2 * bet * U1));
}

double m_12(double U2,double a12) {
    return ((-1) * dt * a12 * ((bet * N - gama) - 2 * bet * U2));
}

double m_21(double U1,double a21) {
    return ((-1) * dt * a21 * ((bet * N - gama) - 2 * bet * U1));
}

double m_22(double U2,double a22) {
    return (1 - dt * a22 * ((bet * N - gama) - 2 * bet * U2));
}


double f(double u){
    return ((bet * N - gama) * u - bet * u * u);
}

void RK2(){
  std::ofstream f1;
  f1.open("RK2.txt");


  double a11=0.25;
  double a12=0.25 - sqrt(3)/6.0;
  double a21=0.25 + sqrt(3)/6.0;
  double a22=0.25;

  double b1=0.5;
  double b2=0.5;

  double c1=0.5 - sqrt(3)/6.0;
  double c2=0.5 + sqrt(3)/6.0;

  double U1=0.,U2=0.;
  double Um1 = 0.0;
  double Um2 = 0.0;
  double dU1=0.,dU2=0.;
  double un = u0;
  double un1=u0;

  for( int i = 0; i <tmax/dt; i++){

    un=un1;
    Um1=0.;
    Um2=0.;

    U1=un;
    U2=un;

    
    
    for(int m=0;m<=iteracje_max;m++){
      if(abs(U1 - Um1)< TOL || abs(U2-Um2)<TOL) break;

      Um1=U1;
      Um2=U2;

      double F1=U1 - un - dt * (a11 * (alfa * U1 - bet * pow(U1,2) + a12 * (alfa * U2 - bet * pow(U2,2))));
      double F2 = U2 - un - dt * (a21 * (alfa * U1 - bet * pow(U1,2) + a22 * (alfa * U2 - bet * pow(U2,2))));

      double m11=m_11(U1,a11);
      double m12=m_12(U2,a12);
      double m21=m_21(U1,a21);
      double m22=m_22(U2,a22);

      dU1=(F2*m12-F2*m22)/(m11*m22-m12*m21);
      dU2=(F2*m21-F2*m11)/(m11*m22-m12*m21);

      U1=Um1+dU1;
      U2=Um2+dU2;
      
    }
    un1=un+dt*(b1*f(U1)+b2*f(U2));
    f1<<i*dt<<" "<<un1<<" "<<N-un1<<"\n";
    
  }
}


int main() {  
  pickardMethod();
  NewtonMethod();
  RK2();
}

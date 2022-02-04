#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

const double delta = 0.01;
const double ro = 1.;
const double mu = 1.;
const int nx = 200;
const int ny = 90;
const int i_1 = 50;
const int j_1 = 55;
const int IT_MAX = 20000;


int check(int i, int j){
    if(i <= i_1 && j <= j_1)
        return 0;
    else
        return 1;
    
}
void WB_psi (double (&psi)[nx+1][ny+1], const double Qwe, const double Qwy, double x[nx+1], double y[ny+1]) {
    // A
    for (int j = j_1; j <= ny; j++) {
        psi[0][j] = (Qwe/(2.*mu)) * ( (pow(y[j], 3)/3.) -  (pow(y[j], 2)/2.)*(y[j_1] + y[ny]) + y[j] * y[j_1] * y[ny] );
    }
    // C
    for (int j = 0; j <= ny; j++) {
        psi[nx][j] = (Qwy/(2.*mu)) * ( (pow(y[j], 3)/3.) -  (pow(y[j], 2)/2.) * y[ny] ) + ( Qwe * pow(y[j_1], 2) * (-y[j_1] + 3. * y[ny]) ) / (12. * mu);
    }
    // B 
    for (int i = 1; i <= nx-1; i++) {
        psi[i][ny] = psi[0][ny];
    }
    //D
    for (int i = i_1; i <= nx-1; i++) {
        psi[i][0] = psi[0][j_1];
    }
    //E    
    for (int j = 1; j <= j_1; j++) {
        psi[i_1][j] = psi[0][j_1];
    }
    //F
    for (int i = 1; i <= i_1; i++) {
        psi[i][j_1] = psi[0][j_1];
    }
}

int border(int i, int j){
    if(i==0 || i == nx)
        return 1;
    if(j==0 || j == ny)
        return 1;
    if(i == i_1 && j <=j_1)
        return 1;
    if(i <= i_1 && j == j_1)
        return 1;
    else
        return 0;
}

void WB_zeta(double (&zeta)[nx+1][ny+1], double psi[nx+1][ny+1], const double Qwe, const double Qwy, double x[nx+1], double y[ny+1]) {
    //A
    for (int j = j_1; j <= ny; j++) {
        zeta[0][j] = (Qwe/(2.*mu)) * (2. * y[j]  - y[j_1] - y[ny]);
    }
    //C
    for (int j = 0; j <= ny; j++) {
        zeta[nx][j] = (Qwy/(2.*mu)) * (2. * y[j] - y[ny]);
    }
    //B
    for (int i = 1; i <= nx-1; i++) {
        zeta[i][ny] = (2./pow(delta, 2)) * (psi[i][ny-1] - psi[i][ny]);
    }
    //D
    for (int i = i_1 + 1; i <= nx-1; i++)  { 
        zeta[i][0] = (2./pow(delta, 2)) * (psi[i][1] - psi[i][0]);
    }
    //E
    for (int j = 1; j <= j_1-1; j++) {
        zeta[i_1][j] = (2./pow(delta, 2)) * (psi[i_1 + 1][j] - psi[i_1][j]);
    }
    //F
    for (int i = 1; i <= i_1; i++) {
        zeta[i][j_1] = (2./pow(delta, 2)) * (psi[i][j_1 + 1] - psi[i][j_1]);
    }
    //wierzcholek E/F
    zeta[i_1][j_1] = (1./2.) * (zeta[i_1 - 1][j_1] + zeta[i_1][j_1 - 1]);
}

double Qwe_count(double Qwe,double y[ny+1]){
  return Qwe * ( ( pow(y[ny], 3) - pow(y[j_1], 3) - 3. * y[j_1] * pow(y[ny], 2) + 3. * pow(y[j_1],2) * y[ny]) / (pow(y[ny],3)) );
}
void solve(const double Qwe) {
      FILE *file1;
  if(Qwe==-1000) file1 = fopen("-1000.dat", "w");
  else if (Qwe==-4000)  file1 = fopen("-4000.dat", "w");
  else  file1 = fopen("4000.dat", "w");
  
  double x[nx+1] = {0.};
    double y[ny+1] = {0.};

    for (int i = 0; i <= nx; i++) {
        x[i]=  delta * i;
    }
    for (int j = 0; j <= ny; j++) {
        y[j] = delta * j;
    }
    double Qwy = Qwe_count(Qwe,y);
    double psi[nx+1][ny+1] = {0.};
    double zeta[nx+1][ny+1] = {0.};
    double u[nx+1][ny+1] = {0.};
    double v[nx+1][ny+1] = {0.};


    WB_psi(psi, Qwe, Qwy, x, y);
    double omega = 0.;
    


    for (int it = 1; it < IT_MAX; it++) {
        if (it < 2000) {
            omega = 0.;
        } else {
            omega = 1.;
        }

        for (int i = 1; i <= nx-1; i++) {
            for (int j = 1; j <= ny-1; j++) {
                if (i > i_1 || j > j_1) {
                    psi[i][j] = (1./4.) * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - pow(delta, 2) * zeta[i][j]);
                    zeta[i][j] = (1./4.) * (zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1]) - omega * (ro/(16.*mu)) * ( (psi[i][j+1] - psi[i][j-1]) * (zeta[i+1][j] - zeta[i-1][j]) - (psi[i+1][j] - psi[i-1][j])*(zeta[i][j+1] - zeta[i][j-1]) );
                }
                u[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2. * delta);
                v[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (2. * delta);
            }
        }



        WB_zeta(zeta, psi, Qwe, Qwy, x, y);

        double Gamma = 0.;
        for (int i = 1; i<= nx-1; i++) {
            Gamma += (psi[i+1][j_1 + 2] + psi[i-1][j_1+2] + psi[i][j_1+3] + psi[i][j_1+1] - 4. * psi[i][j_1+2] - pow(delta,2) * zeta[i][j_1+2]); 
        }
      //  std::cout<<"it="<<it<<"gamma="<<Gamma<<"\n";
      

    } 
    double u1,v1;
    for (int i = 0; i <= nx; i++){
        for (int j = 0; j <= ny; j++) {
          if(!border(i,j)&&check(i,j)){
                u1 = (psi[i][j + 1] - psi[i][j - 1]) / (2 * delta);
                v1 = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * delta);
                 fprintf(file1,"%g %g %g %g %g %g\n",delta*i, delta*j,psi[i][j], zeta[i][j], u1, v1 );
            }
            else{

              fprintf(file1,"%g %g %g %g %g %g\n",delta*i, delta*j,psi[i][j], zeta[i][j], 0., 0. );
            }
    
        }
        fprintf(file1, "\n");

    }   
}

int main() {    
    double Q1 = -1000.;
solve(Q1);
    double Q2 = -4000.;
solve(Q2);
 double Q3 = 4000.;
solve(Q3);
 
    return 0;
}
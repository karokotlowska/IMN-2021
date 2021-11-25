#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;


#define nx  128
#define ny  128

const double TOL = pow(10,-8);
const double delta =0.2;
const double xmax = nx*delta;
const double ymax = ny*delta;

int t[5]={16, 8, 4, 2, 1};


void edge_condition(double V[nx+1][ny+1]){
	for(int y=0; y<ny+1;y++){
	    V[0][y] =(1.0)*sin(M_PI*delta*y/ymax);
	}
    for(int x=0; x<nx+1;x++){
    	V[x][ny]=(-1.0)*sin(2*M_PI*delta*x/xmax);
	}
  	for(int y=0; y<ny+1;y++){
        V[nx][y]=(1.0)*sin(M_PI*delta*y/ymax);
	}
	for(int x=0; x<nx+1;x++){
    	V[x][0]=(1.0)*sin(2*M_PI*delta*x/xmax);
	}
	
}

double is_termination_condition(double V[nx+1][ny+1], int k){
    double S = 0.0;

    for(int i = 0; i <= nx-k; i += k){
		for(int j = 0; j <= ny-k; j += k){
			S += 0.5*pow(k*delta, 2)*(pow((V[i+k][j] - V[i][j])/(2*k*delta) + (V[i+k][j+k]-V[i][j+k])/(2*k*delta), 2) +
					pow((V[i][j+k] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta),2));
		}
	}
    return S;
}

void grid_redefinition(double V[nx+1][ny+1],int k){
  if(k!=1){
	for(int i=0;i<=nx-k;i+=k){
		for(int j=0;j<=ny-k;j+=k){
			V[i+k/2][j+k/2]=0.25*(V[i][j]+V[i+k][j]+V[i][j+k]+V[i+k][j+k]);
			if(i!=nx-k) V[i+k][j+k/2]=0.5*(V[i+k][j]+V[i+k][j+k]);
			if(j!=ny-k) V[i + k/2][j + k] = 0.5*(V[i][j+k] + V[i+k][j+k]);			
			if(j!=0) V[i+k/2][j]=0.5*(V[i][j]+V[i+k][j]);
			if(i!=0) V[i][j+k/2]=0.5*(V[i][j]+V[i][j+k]);    	
	    }
  }
  cout<<k<<" ";
}

}


void solve(ofstream &f1, ofstream &f2){
	double V[nx+1][ny+1];
	for(int i=0;i<nx+1;i++){
		for(int j=0;j<ny+1;j++){
			V[i][j] = 0.0;
		}
	}
	
	edge_condition(V);
	
	double S_prev, S_curr=0.;
	int n=0;
	int sumaiter=0;
  int currentKiter=0;
	for(int k=t[n];n<=(sizeof(t)/sizeof(t[0])),k=t[n];n++){
		do{
      		sumaiter+=1;
          currentKiter+=1;
      		for(int i=k;i<=nx-k; i+=k)
        		for(int j=k;j<=ny-k;j+=k)
					 V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
      		S_prev=S_curr;
      		S_curr=is_termination_condition(V,k);        
      		f1<<" "<<sumaiter<<" "<<S_curr<<"\n";
    	}while(fabs((S_curr-S_prev)/S_prev)>TOL);
    	cout<<k<<" "<<currentKiter<<"\n";
      currentKiter=0;
      f1<<"\n\n";
      for(int i=0;i<nx+1;i+=k){
            for(int j=0;j<ny+1;j+=k){
      		      f2<<delta*i<<" "<<delta*j<<" "<<V[i][j]<<"\n";
            }
        f2<<"\n\n";
      }
        
    	grid_redefinition(V,k);
	}
    
}


int main(){
    ofstream file1;
    file1.open("potential_map.dat");
    ofstream file2;
    file2.open("functional_integral.dat");

    solve(file2, file1);
 

    return 0;
}

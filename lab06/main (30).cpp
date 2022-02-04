#include <algorithm>
#include <cmath>
#include <fstream>
#include "mgmres.h"


#define delta 0.1


int j(int nx,int l){
    return (floor(l/(nx+1)));
}

int i(int nx,int l){
    return l-j(nx,l)*(nx+1);
}

double getEl(int nx,int l,int eps1,int eps2){
   if(i(nx,l) <= nx/2)
        return eps1;
    else
        return eps2;
}


double ro1(double x, double y, double xmax, double ymax, double sig) {
	return exp(-1*pow(x-0.25*xmax,2)/pow(sig,2)-pow(y-0.5*ymax,2)/pow(sig,2));
}

double ro2(double x, double y, double xmax, double ymax, double sig) {
	return (-1)*exp(-1*pow(x-0.75*xmax,2)/pow(sig,2)-pow(y-0.5*ymax,2)/pow(sig,2));
}

double ro(double x, double y, double xmax, double ymax, double sig){
  return ro1(x,y,xmax,ymax,sig)+ro2(x,y,xmax,ymax,sig);
}


int fillSparseMatrix(double V1, double V2, double V3, double V4, int nx, int ny, int *ja, int *ia, double *a, double *b,double xmax, double ymax, int eps1, int eps2, FILE *matrix, FILE *vector, int checkParam) {
	int k=-1;
  int nz_num = 0;
  int N=(nx+1)*(ny+1);
	for(int l=0;l<N;l++){
		int brzeg=0;  
		double vb=0.;  

		if(i(nx, l)==0){
			brzeg=1;
			vb=V1;
		}
    if(j(nx, l)==ny){
			brzeg=1;
			vb=V2;
		}
		if(i(nx,l)==nx){
			brzeg=1;
			vb=V3;
		}
		if(j(nx,l)==0){
			brzeg=1;
			vb=V4;
		}

    b[l]=(-1)*(ro(delta*i(nx, l),delta*j(nx,l),xmax,ymax,xmax/10)); 

		if(brzeg==1){
			b[l]=vb;
    }

		ia[l]=-1; 

		if(l-nx-1>=0&&brzeg==0){
			k++;
			if(ia[l]<0){
        ia[l]=k;
      }
			a[k]=getEl(nx,l,eps1,eps2)/pow(delta,2);
			ja[k]=l-nx-1;
		}

		if(l-1>=0&&brzeg==0) {
			k++;
			if(ia[l]<0){
        ia[l]=k;
      }
			a[k]=getEl(nx,l,eps1,eps2)/pow(delta,2);
			ja[k]=l-1;
		}

		k++;
		if(ia[l]<0){
      ia[l]=k;
    }
		if(brzeg==0){
			a[k]=-(2*getEl(nx,l,eps1,eps2)+getEl(nx,l+1,eps1,eps2)+getEl(nx,l+nx+1,eps1,eps2))/pow(delta,2);
    }
		else{
			a[k]=1;
    }
		ja[k]=l;

		if(l<N&&brzeg==0){
			k++;
			a[k]=getEl(nx,l+1,eps1,eps2)/pow(delta,2);
			ja[k]=l+1;
		}
		if(l<N-nx-1&&brzeg==0){
			k++;
			a[k]=getEl(nx,l+nx+1,eps1,eps2)/pow(delta, 2);
			ja[k]=l+nx+1;
		}
    if(checkParam==0){
      fprintf(vector,"%d %d %d %f \n", l, i(nx,l), j(nx,l), b[l]);
       fprintf(matrix,"%d %d %d %f \n", l, i(nx,l), j(nx,l), a[l]);
      
    }
    }
	  nz_num=k+1; 
	  ia[N]=nz_num;
    return nz_num;
}

void solve(int nx, int ny, int V1, int V2, int V3, int V4,int eps1, int eps2, double xmax, double ymax,  FILE * matrix, FILE * vector,FILE * map,int checkParam){
  int N=(nx+1)*(ny+1);
  double a[5*N]; //niezerowe wartosci elementow macierzowych
	int ja[5*N]; //przechowuje informacje o numerach kolumn
	int ia[N+1]; //wsakzniki do elementow rozpoczynajacych dany wiersz
	double b[N]; //wektor wyrazow wolnych
	double V[N]; //wektor rozwiazan

  for(int i=0;i<N+1;i++){
    ia[i]=-1;
  }

  int nz_num = fillSparseMatrix(V1,V2,V3,V4,nx,ny,ja,ia,a,b,xmax,ymax,eps1,eps2,matrix,vector,checkParam);


  int itr_max=500;
  int mr=500;
  double tol_abs=pow(10,-8);
  double tol_rel=pow(10,-8);

  pmgmres_ilu_cr(N,nz_num,ia,ja,a,V,b,itr_max,mr,tol_abs,tol_rel);


double change=0.;
if(checkParam==1){
	for(int k=0;k<N;k++){
    if(delta*j(nx,k)>change){
			fprintf(map,"\n");
		}
		fprintf(map,"%f %f %f \n",delta*j(nx,k),delta*i(nx, k),V[k]);
		change=delta*j(nx,k);
		}
	}
}


int main(){
	int nx=4;
  int ny=4;

  int V1,V2,V3,V4;
  V1=V3=10;
  V2=V4=-10;

  int eps1,eps2;
  eps1=eps2=1;

  double xmax=0;
  double ymax=0;

 	FILE *matrix = fopen("matrixCheck.txt", "w");
  FILE *vector = fopen("vectorCheck.txt", "w");	
  FILE *map = fopen("mapEmpty.txt", "w");
  solve(nx, ny, V1, V2, V3, V4,eps1, eps2, xmax, ymax,  matrix, vector, map,0);
  fclose(matrix);
  fclose(vector);
  fclose(map);

  matrix=fopen("txt.txt","w");
  vector=fopen("txt.txt","w");
  map=fopen("map_n_50.txt","w");
  nx=ny=50;
  solve(nx, ny, V1, V2, V3, V4,eps1, eps2, xmax, ymax, matrix, vector, map,1);
  fclose(map);


	nx=ny=100;
	map = fopen("map_n_100.txt", "w");
  solve(nx, ny, V1, V2, V3, V4,eps1, eps2, xmax, ymax,  matrix,vector, map,1);
  fclose(map);
	

	nx=ny=200;	
	map = fopen("map_n_200.txt", "w");
    solve(nx, ny, V1, V2, V3, V4,eps1, eps2, xmax, ymax,  matrix, vector, map,1);
    fclose(map);



//---------------------------------------------------------


	nx=ny=100;
	V1=V2=V3=V4=0;
	xmax=delta*nx;
	ymax=delta*ny;
  
 eps1=eps2=1;
  vector=fopen("txt.txt","w");
	map = fopen("map_eps_1.txt", "w");
  solve(nx, ny, V1, V2, V3, V4,eps1, eps2, xmax, ymax,  matrix, vector, map,1);
  fclose(map);

	eps1 = 1;
	eps2 = 2;
	map = fopen("map_eps_1_2.txt", "w");
  solve(nx, ny, V1, V2, V3, V4,eps1, eps2, xmax, ymax,  matrix,vector, map,1);
  fclose(map);


	eps1 = 1;
	eps2 = 10;
	map = fopen("map_eps_1_10.txt", "w");
  solve(nx, ny, V1, V2, V3, V4,eps1, eps2, xmax, ymax,  matrix, vector, map,1);
  fclose(map);
  fclose(matrix);
  fclose(vector);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//zastawka miedzy 200 210 (po x) i 50 (po y)
#define j1 50
int nx = 400;
int ny = 90;
int i1 = 200;
int i2 = 210;
double delta = 0.01;
double sigma = 0.1;
double xA = 0.45;
double yA = 0.45;
int IT_MAX = 10000; 




int zastawka(int i, int j){
	if(i >= i1 && i <= i2 && j<= j1)
		return 1;
	return 0;
}

void solve(double D, double** v_x, double** v_y, double dt, FILE** f1, FILE** f2){
	double** u0 = calloc(nx+1, sizeof(double*));
	for(int i = 0; i <= nx; i++)	{
		u0[i] = calloc(ny+1, sizeof(double));
	}

	double** u1 = calloc(nx+1, sizeof(double*));
	for(int i = 0; i <= nx; i++)	{
		u1[i] = calloc(ny+1, sizeof(double));
	}

	//inicjalizacja gestosci u0 - wzor 17 //warunek pocztkowy //Xa Ya to polozenie srodka gestosci a sigma to rozmycie
	for(int i = 0; i <= nx; i++){
		for(int j = 0; j <= ny; j++){
			u0[i][j] = exp( -(pow(i*delta - xA, 2) + pow(j*delta - yA, 2))/(2.0*pow(sigma, 2)) )
				/(2.0*3.14*pow(sigma, 2));
		}
	}

	for(int IT = 1; IT <= IT_MAX; IT++)	{
		//inicjalizacja kolejnego kroku u1 = u0
		for(int i = 0; i <= nx; i++)		{
			for(int j = 0; j <= ny; j++)			{
				u1[i][j] = u0[i][j];
			}
		}
		
//algorytm dla rwonania AD
		for(int k = 1; k <= 20; k++){
			for(int i = 0; i <= nx; i++){
				for(int j = 1; j <= ny-1; j++){					
					if(zastawka(i,j))
					{
						continue;
					}
					else if(i == 0 || i == nx)
					{
						if(i == 0)
						{
							u1[i][j] = (1.0/(1.0 + 2.0*D*dt/pow(delta, 2)))
								*(u0[i][j] - 0.5*dt*v_x[i][j]*
								(0.5*(u0[i+1][j] - u0[nx][j])/delta + (0.5*(u1[i+1][j] - u1[nx][j])/delta) ) 
								- 0.5*dt*v_y[i][j] * ( (0.5*(u0[i][j+1] - u0[i][j-1])/delta) + (0.5*(u1[i][j+1] - u1[i][j-1])/delta) ) 
								+ 0.5*D*dt * ( (u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta, 2)
								+ (u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1])/pow(delta, 2) ) );
						}
						if(i == nx)
						{
							u1[i][j] = (1.0/(1.0 + 2.0*D*dt/pow(delta, 2)))
								*(u0[i][j] - 0.5*dt*v_x[i][j]*
								(0.5*(u0[0][j] - u0[i-1][j])/delta + (0.5*(u1[0][j] - u1[i-1][j])/delta) ) 
								- 0.5*dt*v_y[i][j] * ( (0.5*(u0[i][j+1] - u0[i][j-1])/delta) + (0.5*(u1[i][j+1] - u1[i][j-1])/delta) ) 
								+ 0.5*D*dt * ( (u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta, 2)
								+ (u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1])/pow(delta, 2) ) );
						}
					}
					else
					{
						u1[i][j] = (1.0/(1.0 + 2.0*D*dt/pow(delta, 2)))
							*(u0[i][j] - 0.5*dt*v_x[i][j]*
							(0.5*(u0[i+1][j] - u0[i-1][j])/delta + (0.5*(u1[i+1][j] - u1[i-1][j])/delta) ) 
							- 0.5*dt*v_y[i][j] * ( (0.5*(u0[i][j+1] - u0[i][j-1])/delta) + (0.5*(u1[i][j+1] - u1[i][j-1])/delta) ) 
							+ 0.5*D*dt * ( (u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta, 2)
							+ (u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1])/pow(delta, 2) ) );
					}
				}
			}
		}

		for(int i = 0; i <= nx; i++){
			for(int j = 0; j <= ny; j++){
				u0[i][j] = u1[i][j];   //zachowuje rozwiazania
			}
		}

		//wyznaczanie c i x_sr 
		double c = 0.;
		double x_sr = 0.;
		for(int i = 0; i <= nx; i++){
			for(int j = 0; j <= ny; j++){
				c += u0[i][j]; ///wzor 19
				x_sr += delta*i * u0[i][j]; //wzor 20
      }
    }

    //wpisywanie
   fprintf(*f2,  "%g %g\n",c*pow(delta, 2), x_sr*pow(delta, 2));
for(int i=0;i<=nx;i++){
  for(int j=0;j<=ny;j++){
    if(IT%(IT_MAX/5)==0){
      fprintf(*f1, "%g %g %g\n", i*delta, j*delta, u0[i][j]);
    }
    else{
      fprintf(*f1,"\n\n");
    }
  }
}
	for(int i = 0; i <= nx; i++)
	{	
		free(u0[i]);
		free(u1[i]);
	}
	free(u0);
	free(u1);

  } 
}


int main()
{
	FILE* f_psi = fopen("psi.dat", "r");
	FILE* fp_u0 = fopen("D0.txt", "w+");
	FILE* fp_u1 = fopen("D1.txt", "w+");
	FILE* fp1 = fopen("c_v_x.txt", "w+");
	
	double **psi=calloc(nx+1,sizeof(double));
  double** v_x = calloc(nx+1, sizeof(double*));
 	double** v_y = calloc(nx+1, sizeof(double*));

	for(int i = 0; i <= nx; i++){
		psi[i] = calloc(ny+1, sizeof(double));
	}
	for(int i = 0; i <= nx; i++){
		v_x[i] = calloc(ny+1, sizeof(double));
	}

	for(int i = 0; i <= nx; i++){
		v_y[i] = calloc(ny+1, sizeof(double));
	}

	//wczytanie psi z pliku
	for(int i = 0; i <= nx; i++){
		for(int j = 0; j <= ny; j++){
			long int a, b;
			fscanf(f_psi, "%ld %ld %lf", &a, &b, &psi[i][j]);
		}
	}

	//wyznaczenie pol predkosci 2.4
	for(int i = 1; i <= nx-1; i++){
		for(int j = 1; j <= ny-1; j++){
			//zastąpienie pochodnych symetrycznymi ilorazami róznicowymi - wzór 11 i 12
				v_x[i][j] = 0.5*(psi[i][j+1] - psi[i][j-1])/delta;
				v_y[i][j] = -0.5*(psi[i+1][j] - psi[i-1][j])/delta;
			}
		}

for(int i = 1; i <= nx-1; i++){
		for(int j = 1; j <= ny-1; j++){
			if(i >= i1 && i <= i2 && j <= j1)
			{
				v_x[i][j] = v_y[i][j] = 0.0; //ustawienie na zastawce wzor 13
			}
		}
    //ustawienie na obu brzegach, dolnym i górnym
		v_x[i][0]=v_y[i][ny] = 0.0;
	}

	for(int j = 0; j <= ny; j++)  //prawy i lewy brzeg - przepisuje predkości z sąsiednich węzłów
	{
		v_x[0][j] = v_x[1][j];
		v_x[nx][j] = v_x[nx-1][j]; //wzor 15
		v_y[0][j] = v_y[1][j];
		v_y[nx][j] = v_y[nx-1][j]; //wzor 16
	}

	FILE* fp_v = fopen("v.txt", "w+");
	for(int i = 0; i <= nx; i++){
		for(int j = 0; j <= ny; j++){
			fprintf(fp_v, "%g %g %g %g\n", i*delta, j*delta, v_x[i][j], v_y[i][j]);
		}
		fprintf(fp_v, "\n");
	}
	fclose(fp_v);

	//wyznaczanie vmax i dt - szukam punktu i,j dla ktorego modul predkosci jest najwiekszy
  double temp = sqrt(pow(v_x[0][0], 2) + pow(v_y[0][0], 2));
	for(int i = 0; i <= nx; i++){
		for(int j = 0; j <= ny; j++){
			double val = sqrt(pow(v_x[i][j], 2) + pow(v_y[i][j], 2));
			if(temp <= val)
				temp = val;
		}
	}

	double v_max = temp;
	double dt = 0.25*delta/v_max; //wzor 18
	
	//algorytm adwekcji dyfuzji
	solve(0.0, v_x, v_y, dt, &fp_u0, &fp1);
	solve(0.1, v_x, v_y, dt, &fp_u1, &fp1);

	fclose(fp_u0);
	fclose(fp_u1);
	fclose(fp1);

	for(int i = 0; i <= nx; i++)
	{
		free(psi[i]);
		free(v_x[i]);
		free(v_y[i]);
	}
	free(psi);
	free(v_x);
	free(v_y);
	
	return 0;
}
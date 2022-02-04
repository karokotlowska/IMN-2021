#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>



#define n_x 40
#define n_y 40
#define N ((n_x+1)*(n_y+1))

#define T_A 40
#define T_B 0
#define T_C 30
#define T_D 0

#define k_B 0.1
#define k_D 0.6

#define delta 1
#define delta_t 1

#define IT_MAX 2000


double pow_T(gsl_vector* T, int l ){
	return ( gsl_vector_get(T,l+1) -2 * gsl_vector_get(T,l) + gsl_vector_get(T,l-1) ) / pow( delta, 2. ) + ( gsl_vector_get(T,l+n_x+1) -2 * gsl_vector_get(T,l) + gsl_vector_get(T,l-n_x-1) ) / pow( delta, 2. );
}


int get_l(int i, int j) {
    return i + j * (n_x + 1);
}


int get_i(int j, int l) {
    return l - j * (n_x + 1);
}

int get_j(int l) {
    return l / (n_x + 1);
}

void solve(){
	FILE *file_T = fopen( "T.dat", "w");
	FILE *file_DT = fopen( "DT.dat", "w");
	
 	int l;
	int signum = 0;
	
	
	gsl_matrix *A = gsl_matrix_calloc(N, N);
	gsl_matrix *B = gsl_matrix_calloc(N, N);
	gsl_vector *c= gsl_vector_calloc(N);
        gsl_permutation *p=gsl_permutation_alloc(N);
	gsl_vector *d = gsl_vector_calloc(N);
	gsl_vector *T= gsl_vector_calloc(N);
	

	for(int i = 1; i < n_x; i++){
		for(int j = 1; j < n_y; j++){
			
			l = get_l(i,j);
			
			gsl_matrix_set( A, l, l-n_x-1, delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( A, l, l-1,     delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( A, l, l+1,     delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( A, l, l+n_x+1, delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( A, l, l,      -(2 * delta_t) / ( pow( delta, 2) ) - 1 );
			
			gsl_matrix_set( B, l, l-n_x-1, -delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( B, l, l-1,     -delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( B, l, l+1,     -delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( B, l, l+n_x+1, -delta_t / ( 2 * pow( delta, 2) ) );
			gsl_matrix_set( B, l, l,        (2 * delta_t) / ( pow( delta, 2) ) - 1 );

		}
	}

	// prawy i lewy brzeg
	int i = 0;
	for(int j = 0; j <= n_y; j++){
		l = get_l(i,j);
		
		gsl_matrix_set( A, l, l, 1 );
		gsl_matrix_set( B, l, l, 1 );
		gsl_vector_set(c, l, 0);
	}
	
	i = n_x;
	for(int j = 0; j <= n_y; j++){
		l = get_l(i,j);
		
		gsl_matrix_set( A, l, l, 1 );
		gsl_matrix_set( B, l, l, 1 );
		gsl_vector_set(c, l, 0);
	}
	
	//gorny brzeg
	int j = n_y;
	for(int i = 1; i < n_x; i++){
		
		l = get_l(i,j);
		
		gsl_matrix_set( A, l, l-n_x-1,  - 1 / (k_B * delta) );
		gsl_matrix_set( A, l, l,        1 + 1 / (k_B * delta) );
		gsl_vector_set( c, l, T_B);
	}
	
	for( int k = 0; k < N; k++ ){
		for(int i = 1; i < n_x; i++){
		
			l = get_l(i,j);
			
			gsl_matrix_set( B, l, k, 0 );
		}	
	}
	
	// dolny rzeg
	j = 0;
	for(int i = 1; i < n_x; i++)
	{
		l = get_l(i,j);
		
		gsl_matrix_set( A, l, l,       1 + 1 / (k_D * delta) );
		gsl_matrix_set( A, l, l+n_x+1, - 1 / (k_D * delta) );
		gsl_vector_set( c, l, T_D);
		
	}
	for( int k = 0; k < N; k++ ){
		for(int i = 1; i < n_x; i++){
		
			l = get_l(i,j);
			
			gsl_matrix_set( B, l, k, 0 );
		}	
	}
	
	
	i = 0;
	for(int j = 0; j < n_y; j++)	{
		l = get_l(i,j);
		gsl_vector_set( T, l, T_A);
	}
	i = n_x;
	for(int j = 0; j < n_y; j++){
		l = i + j * ( n_x + 1 );
		gsl_vector_set( T, l, T_C);
	}
	for(i = 1; i <= n_x - 1; i++){
        for(int j = 0; j <= n_y; j++){
			l = get_l(i,j);
			gsl_vector_set(T, l, 0);
		}
  	}
	
	gsl_linalg_LU_decomp ( A , p , &signum );
	
	for( int it = 1; it <= IT_MAX; it++){
		gsl_blas_dgemv( CblasNoTrans , 1, B, T, 0, d );
		gsl_blas_daxpy ( 1 , c , d );
		gsl_linalg_LU_solve ( A , p, d , T );
		
		
		if( it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000 ){
    		for(int i = 1; i < n_x; i++){
				  for(int j = 1; j < n_y; j++){
					l = get_l(i,j);
					fprintf(file_T, "%d\t %d\t %d\t %g \n", it, i * delta, j * delta, gsl_vector_get(T, l) );
					fprintf(file_DT, "%d\t %d\t %d\t %g \n", it, i * delta, j * delta, pow_T(T,l) );
				}
				fprintf(file_T, "\n");
				fprintf(file_DT, "\n");
			}
			fprintf(file_T, "\n");
			fprintf(file_DT, "\n");
		}		
		
	}

	gsl_matrix_free(A);
 	gsl_vector_free(T);
	gsl_vector_free(d);
	gsl_permutation_free(p);
	gsl_matrix_free(B);
	gsl_vector_free(c);
	fclose(file_T);
	fclose(file_DT);
}


int main()
{
	solve();
}



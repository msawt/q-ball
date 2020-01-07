#include "DirOfInterest.h"
#include "Matrices.h"
#include <math.h>
#include<stdlib.h>

#define pi 3.14159265

double* getEquator(const int k){
	double **C = (double** )malloc(sizeof(double *)*k);
	for(int i=0; i < k; i++) {
	  C[i] = (double *)malloc(sizeof(double)*3);
	}
	for(int i = 0; i < k; ++i){
		C[0][i] = cos((2 * pi / k) * (i + 1));
		C[1][i] = sin((2 * pi / k) * (i + 1));
		C[2][i] = 0;
	}
	return *C;
}

double* getRotationMat(double z[3],double u[3]){
	double **prod = (double** )malloc(sizeof(double *)*k);
		for(int i=0; i < k; i++) {
			prod[i] = (double *)malloc(sizeof(double)*3);
		}

	double add[3];

	add[0] = z[0] + u[0];
	add[1] = z[1] + u[1];
	add[2] = z[2] + u[2];

	matMulPoints(add, add, prod);
	double div = mulPoints(z, u, 3) + 1;

	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			prod[i][j] /= div;
		}
	}
	prod[0][0] --; //Subtract I from the matrix
	prod[1][1] --;
	prod[2][2] --;

	return * prod;
}

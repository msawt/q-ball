#include "DirOfInterest.h"
#include "Matrices.h"
#include <math.h>
#include<stdlib.h>
#include<stdio.h>

#define pi 3.14159265

void getEquator(const int k, double ** C){

	//In:
	//	k sampling directions
	//	3 x k output matrix C

	//Operation:
	//	Takes k sampling directions and returns an equator of k equally spaced points

	//C must be 3 x k
	//To allocate space for c (use in main):

	for(int i = 0; i < k; ++i){
		C[0][i] = cos((2 * pi / k) * (i + 1));
		C[1][i] = sin((2 * pi / k) * (i + 1));
		C[2][i] = 0;
	}
}

void getRotationMat(double* z,double* u, double** out){

	//In:
	//	3 x 1 corresponding equator point z
	//	3 x 1 direction of interest u
	//	3 x 3 output rotation matrix for that direction

	//Operation:
	//	Returns the matrix that you multiply each equator point by in order to rotate the equator to u



	double add[3];
	double prod[3][3];

	add[0] = z[0] + u[0];
	add[1] = z[1] + u[1];
	add[2] = z[2] + u[2];

	matMulPoints(add, add, prod);
	double div = mulPoints(z, u, 3) + 1;

	//printf("Denom: %f\n", div);

	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			prod[i][j] = prod[i][j] / div;
		}
	}
	prod[0][0] --; //Subtract I from the matrix
	prod[1][1] --;
	prod[2][2] --;

	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++j){
			out[i][j] = prod[i][j];
		}
	}

}

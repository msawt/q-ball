#include"ConvolutionAndPrediction.h"
#include"Interpolation.h"
#include<stdlib.h>
#include<stdio.h>
#include"Matrices.h"

void getDiffusionSignal(double **q, double **v, int qLen, int vLen,double**out){
	//INPUT: Two arrays of points, and their respective lengths
	//Out: Two dimensional m x p array H, dynamically allocated

	//NOTE: Also used in finding diffusion signal estimate

	//printf("getDiffusionSignal\n");

	for(int i = 0; i < qLen; ++i){
		double* tempQ = malloc(sizeof(double)*3);
		tempQ[0] = q[0][i];
		tempQ[1] = q[1][i];
		tempQ[2] = q[2][i];

		//printf("Q[%d] = %f, %f, %f\n",i,tempQ[0],tempQ[1],tempQ[2]);

		for(int j = 0; j < vLen; ++j){

			//printf("i = %d, j = %d\n", i, j)


			double* tempV = malloc(sizeof(double)*3);
			tempV[0] = v[0][j]; 
			tempV[1] = v[1][j]; 
			tempV[2] = v[2][j];

			double z = getKernel(tempQ,tempV);

			out[i][j] = z;


			free(tempV);
			//printf("\n");
		}
		free(tempQ);
	}

}

double* estimateWeightVector(double ** H, int rows, int cols, double * e, int eLen) {
	//Dimensions of H: m x p
	//INPUT: q x v matrix H, vector e with length eLen
	//OUTPUT: Array w, estimated weight vector
	//Multiplies the pseudoinverse of H by e.

	double ** HPlus = malloc(sizeof(double*) * rows);
	for(int i = 0; i < rows; ++i){
			HPlus[i] = malloc(sizeof(double) * cols);
		}

	pseudoInverse(H,rows, cols, HPlus);

	//HPlus Dimensions: cols x rows (p x m)
	//e Dimensions: m x 1

	double* out = malloc(sizeof(double)*cols);
	for(int i = 0; i < cols; ++i){ //Go through every ROW in HPlus
		double sum = 0;
		for(int j = 0; j < rows; ++j) { //Go through every COL in HPlus and multiply it by e[j]
			sum += HPlus[j][i] * e[j];
		}
		out[i] = sum;
	}
	return out;
}

void getReconstructionMatrix(double ** G, double ** H, int k, int n, int p, int m, double ** out) {
	//Result: (In X 1kT)GH+
		//X is the matrix direct product
	//the summation should be implemented by repartitioning the (kn) x m matrix GH+ into a k x n x m array and then summing over the first dimension.

	//Dimensions of G: (kn) x p
	//Dimensions of H+: p x m
	//Dimensions of GH+: kn x m

	//Dimensions of out:
		//n x m

	double ** mul = malloc(sizeof(double*) * k * n);
	for(int i = 0; i < k * n; ++ i){
		mul[i] = malloc(sizeof(double) * m);
	}

	//Calculate GH+
	for(int row = 0; row < (k*n); row++){
		for(int col = 0; col < m; col++){
			mul[row][col] = 0;
			for(int mark = 0; mark < p; mark++){
				mul[row][col] += G[row][mark] * H[mark][col];
			}
		}
	}

	/*
	printf("Mul: \n");
	for(int i = 0; i < (k*n);++i){
		printf("%lf\t%lf\t%lf\t%lf\n",mul[i][0],mul[i][1],mul[i][2],mul[i][3]);
	}
	*/

	//n = 2
	//k = 3
	//Then sum over equator
	//k*n rows to k rows
	for(int col = 0; col < m; col++){ //For every column,
		for(int row = 0; row < n*k; row+=k){ //Iterate through every row and add k rows into 1 row
			double sum = 0;
			for(int j = 0; j < k; j++){ //Add the row values from i to k into one figure, and then set the output to it
				sum += mul[row+j][col];
			}
			out[row/k][col] = sum;
		}
	}

	for(int i = 0; i < (k * n); ++ i){
		free(mul[i]);
	}
	free(mul);

}


double calculateNormalizationConstant(double ** A, double * e, int n, int m){
	//Calculates the normalization constant Z = 1nTAe

	// 1nT = 1 x n vector of ones
	//A = Reconstruction matrix (use getReconstructionMatrix)
		//n x m
	//e = Diffusion Signal


	//Calculate 1nT*A
	//Resulting vector: 1 x m
	//(Sum each col)

	//printf("calculateNormalizationConstant: \n");

	double firstProduct[m];
	for(int i = 0; i < m; ++i){
		double sum = 0;
		for(int j = 0; j < n; ++j){
			sum += A[j][i];
		}
		firstProduct[i] = sum;
	}

	/*
	printf("First Matrix Product:\n");
	for(int i = 0; i < m; ++i) {
		printf("%lf\t",firstProduct[i]);
	}
	printf("\n");

	*/

	//Calculate firstProduct * e
	double out = 0;
	for(int i = 0; i < m; ++i) {
		out += e[i] * firstProduct[i];
	}

	return out;
}

void computeODF(double ** A, double * e, int n, int m, double * out) { //Normalization constant might not be necessary, since we are looking for peaks rather than values
	//Inputs: Reconstruction matrix A, diffusion signal e

	//Calculates the estimated ODF, (1/Z)Ae, with Z being a normalization constant

	//Output: (1/Z)Ae, n x 1 vector


	//Normalization Constant
	double norm = 1 / calculateNormalizationConstant(A,e,n,m);

	//Calculate Ae and multiply each term by norm
	for(int i = 0; i < n; ++i){
		double count = 0;
		for(int j = 0; j < m; ++j){
			count += e[j] * A[i][j];
		}
		out[i] = count*norm;
		//printf("ODF[%d]: %f\t",i,count);
	}
	//printf("\n");

}



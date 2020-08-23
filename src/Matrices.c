#include"Matrices.h"
#include<stdlib.h>
#include<math.h>
#include<stdio.h>

double mulPoints(double* rowM, double* colM, int length){
	//printf("mulPoints\n");
	double out = 0;
	for(int i = 0; i < length; ++i){
		out += rowM[i] * colM[i];
		//printf("\t\tmulPoints i = %d\n", i);
	}
	if(isnan(out)){
		printf("Expression is nan. Inputs:\t");
		for(int i = 0; i < length; ++i){
			printf("%f\t",rowM[i]);
		}
		printf(" and ");
		for(int i = 0; i < length; ++i){
			printf("%f\t",colM[i]);
		}
		printf("\n");
	}
	return out;
}

void matMulPoints(double p1[3], double p2[3], double out[3][3]){

	out[0][0] = p1[0] * p2[0];
	out[0][1] = p1[0] * p2[1];
	out[0][2] = p1[0] * p2[2];

	out[1][0] = p1[1] * p2[0];
	out[1][1] = p1[1] * p2[1];
	out[1][2] = p1[1] * p2[2];

	out[2][0] = p1[2] * p2[0];
	out[2][1] = p1[2] * p2[1];
	out[2][2] = p1[2] * p2[2];
}

void div33ByScalar(double m1[3][3], double div, double out[3][3]){
	out[0][0] = m1[0][0] / div;
	out[0][1] = m1[0][1] / div;
	out[0][2] = m1[0][2] / div;

	out[1][0] = m1[1][0] / div;
	out[1][1] = m1[1][1] / div;
	out[1][2] = m1[1][2] / div;

	out[2][0] = m1[2][0] / div;
	out[2][1] = m1[2][1] / div;
	out[2][2] = m1[2][2] / div;
}

void mulThreeByThree(double m1[3][3], double m2[3][3], double out[3][3]){
	out[0][0] = (m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0]);
	out[0][1] = (m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1] + m1[0][2] * m2[2][1]);
	out[0][2] = (m1[0][0] * m2[0][2] + m1[0][1] * m2[1][2] + m1[0][2] * m2[2][2]);

	out[1][0] = (m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0] + m1[1][2] * m2[2][0]);
	out[1][1] = (m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1] + m1[1][2] * m2[2][1]);
	out[1][2] = (m1[1][0] * m2[0][2] + m1[1][1] * m2[1][2] + m1[1][2] * m2[2][2]);

	out[2][0] = (m1[2][0] * m2[0][0] + m1[2][1] * m2[1][0] + m1[2][2] * m2[2][0]);
	out[2][1] = (m1[2][0] * m2[0][1] + m1[2][1] * m2[1][1] + m1[2][2] * m2[2][1]);
	out[2][2] = (m1[2][0] * m2[0][2] + m1[2][1] * m2[1][2] + m1[2][2] * m2[2][2]);
}

double det(double **a,int n)
{
   int i,j,j1,j2;
   double d = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      d = a[0][0];
   } else if (n == 2) {
      d = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      d = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++) {
            m[i] = malloc((n-1)*sizeof(double));
         }
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         d += pow(-1.0,1.0+j1+1.0) * a[0][j1] * det(m,n-1);
         for (i=0;i<n-1;i++){
            free(m[i]);
         }
         free(m);
      }
   }
   /*
   printf("Det of \n");
   for(int row = 0; row < n; ++row){
	   for(int col = 0; col < n; ++col){
		   printf("%lf   ",a[row][col]);
	   }
	   printf("\n");
   }
   printf("n: %d \n",n);
   printf("Det: %lf \n",d);
	*/
   return(d);
}

void getCfactor(double ** M, double ** t,int n) {
   if(n==2){
	   t[0][0] = M[1][1];
	   t[0][1] = M[1][0];
	   t[1][0] = M[0][1];
	   t[1][1] = M[0][0];
   }
   else{
	   for(int i = 0; i < n; ++i){
		   for(int j = 0; j < n; ++j){
			   double ** submatrix = malloc(sizeof(double*)* n-1);
			   for(int _ = 0; _ < n-1; ++_){
				   submatrix[_] = malloc(sizeof(double)*n-1);
			   }
			   int p = 0;
			   int q = 0;
			   for(int x = 0; x < n; ++x){
				   for(int y = 0; y < n; ++y){
					   if(x != i && y != j){
						   submatrix[p][q] = M[i][j];
						   if(q < n-2){
							   ++q;
						   }
						   else{
							   q = 0;
							   ++p;
						   }
					   }
				   }
			   }

			   t[i][j] = det(submatrix,n-1);

			   for(int _ = 0; _ < n-1; ++_){
				   free(submatrix[_]);
			   }
			   free(submatrix);
		   }
	   }
   }
   /*
   printf("Original: \n");
   for(int row = 0; row < n; ++row){
   	   for(int col = 0; col < n; ++col){
   		   printf("%lf   ",M[row][col]);
   	   }
   	   printf("\n");
      }
   printf("Cofactor: \n");
      for(int row = 0; row < n; ++row){
      	   for(int col = 0; col < n; ++col){
      		   printf("%lf   ",t[row][col]);
      	   }
      	   printf("\n");
         }
         */
}

void ADJ(double ** M,double ** adj, int N) {
   if (N == 1) {
      adj[0][0] = 1; return;
   }
   //Allocate space for cofactor matrix

   double ** t = malloc(sizeof(double*)*N);
   for(int i = 0; i < N; ++i){
	   t[i] = malloc(sizeof(double) * N);
   }

   getCfactor(M,t,N);

   for(int row = 0; row < N; ++row){
   	   for(int col = 0; col < N; ++col){
   		   if((row+col+2) % 2 == 0){
   			   adj[row][col] = t[col][row];
   		   }
   		   else{
   			   adj[row][col] = -1 * t[col][row];
   		   }
   	   }
      }

   for(int i = 0; i < N;++i){
	   free(t[i]);
   }
   free(t);
}


void pseudoInverse(double ** H, int rowDim, int colDim, double ** out) {
	//Computes the Moore-Penrose pseudoinverse of H, H+

	//H+ = (Ht * H)^-1 * Ht

	//Input: m x n
	//Ht * H: n x n
	//Inverse: n x n
	//Output: n x m

	//Define H transpose

	double Ht[colDim][rowDim];

	for(int i = 0; i < rowDim; ++i){
		for(int j = 0; j < colDim; ++j){
			Ht[j][i] = H[i][j];
		}
	}

	double** mul1 = malloc(sizeof(double*) * colDim);
	    for(int i = 0; i < colDim; ++i){
	    	mul1[i] = malloc(sizeof(double) * colDim);
	    }

	printf("Ht: \n");
		for(int row = 0; row < colDim; ++row){
			for(int col = 0; col < rowDim; ++col){
				printf("%lf   ",Ht[row][col]);
			}
		printf("\n");
		}
	printf("\n");

	//Compute Ht * H
	//Ht: colDim x rowDim H: rowDim x colDim

	for(int row = 0; row < colDim; row++){
		for(int col = 0; col < colDim; col++){
			mul1[row][col] = 0;
			for(int mark = 0; mark < rowDim; mark++){
				mul1[row][col] += Ht[row][mark] * H[mark][col];
			}
		}
	}

	printf("Ht * H: \n");
		for(int row = 0; row < colDim; ++row){
			for(int col = 0; col < colDim; ++col){
				printf("%lf   ",mul1[row][col]);
			}
		printf("\n");
		}
	printf("\n");


	//Compute Determinant
	double determinant = det(mul1,colDim);

	printf("Determinant of Ht * H: %lf\n", determinant);

	//Allocate space for adjoint
	double** adjoint = malloc(sizeof(double*) * colDim);
		for(int i = 0; i < colDim; ++i){
			adjoint[i] = malloc(sizeof(double) * colDim);
		}

	ADJ(mul1,adjoint,colDim);

	printf("Adjoint: \n");
		for(int row = 0; row < colDim; ++row){
			for(int col = 0; col < colDim; ++col){
				printf("%lf   ",adjoint[row][col]);
			}
		printf("\n");
		}
	printf("\n");

	//Free up space for mul1
	for(int i = 0; i < colDim;++i){
		free(mul1[i]);
	}
	free(mul1);

	//Inverse = 1 / det * adjoint

	double invDet = 1 / determinant;
	for(int i = 0; i < colDim; ++i){
		for(int j = 0; j < colDim; ++j){
			adjoint[i][j] *= invDet;
		}
	}


	//Multiply by Ht
	for(int row = 0; row < colDim; row++){
		for(int col = 0; col < rowDim; col++){
			out[row][col] = 0;
			for(int mark = 0; mark < colDim; mark++){
				out[row][col] += adjoint[row][mark] * H[col][mark];
			}
		}
	}

	//Free up adjoint
	for(int i = 0; i < colDim;++i){
		free(adjoint[i]);
	}
	free(adjoint);
}




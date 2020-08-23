#include"Interpolation.h"
#include"Matrices.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double getKernel(double* n1, double* n2){
	//printf("getKernel\n");
	int width = 1; //WIDTH PARAMETER, LEFT CONSTANT FOR NOW

	//printf("mulPoints: %f\n", mulPoints(n1, n2, 3));
	double d = acos(abs(mulPoints(n1, n2, 3)));
	//printf("d: %f\n", d);
	//printf("Finished calling mulPoints\n");
	double out = exp(-pow(d,2)/pow(width,2));
	//printf("kernel: %f\n",out);
	return out;
}

#include"Interpolation.h"
#include"Matrices.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double getKernel(double* n1, double* n2){
	//printf("getKernel\n");
	double width = 1; //WIDTH PARAMETER, LEFT CONSTANT FOR NOW

	double mul = mulPoints(n1,n2,3);

	if(mul < 0){
		mul *= -1;
	}

	if(mul > 1 && mul <= 1.0001){ //Catch any floating point errors
		mul = 1;
	}
	else if(mul < -1 && mul >= -1.0001){
		mul = -1;
	}

	double d = acos(mul);


	//spherical Gaussian function
	
	double nd2 = (d*d)*-1;
	double out = exp(nd2/(width*width));
	
	

	//inverse multiquadric function
	//double out = 1/sqrt((d*d)+(width*width));

	//printf("\t\tkernel: %f\n",out);
	return out;
}

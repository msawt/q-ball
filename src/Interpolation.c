#include"Interpolation.h"
#include<math.h>

double getKernel(double n1[3], double n2[3]){
	int width = 1; //WIDTH PARAMETER, LEFT CONSTANT FOR NOW
	double d = acos(mulPoints(*n1, *n2, 3));
	return exp(-d^2/width^2);
}

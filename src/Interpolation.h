/*
 * First, specify a distance metric d(n1,n2), which is equal to the minimum angle between
 * direction vectors.
 * Next, the spherical Gaussian interpolation kernel is implemented, though there are
 * others that can be used.
 * 		- Thing(a) = exp(-a^2/s^2), where a = d(n1,n2) and s is a width parameter.
 * 		- Width controls tradeoff between accuracy and stability.
 * 		- Small width = high accuracy, vice versa for large width
 * 		- Possible to find optimal value (Check later in the paper)
 *
 *
 *
 */


#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

double getKernel(double* n1, double* n2);



#endif /* INTERPOLATION_H_ */

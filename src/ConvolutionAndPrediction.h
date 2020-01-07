/*
Then, specify p unit vectors {v} (sampling directions). The diffusion signal
 * 	can be expressed as a convolution of the spherical basis functions
 * 		- e = Hw
 * 		- w is the coefficient vector
 * 		- H(ij) = thing(d(qi, vj)
 *
*/

double* getDiffusionSignal(double *q, double *v, int qLen, int vLen);

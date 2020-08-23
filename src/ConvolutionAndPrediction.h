/*
Then, specify p unit vectors {v} (sampling directions). The diffusion signal
 * 	can be expressed as a convolution of the spherical basis functions
 * 		- e = Hw
 * 		- w is the coefficient vector
 * 		- H(ij) = X(d(qi, vj)
 *
*/

#ifndef ConvolutionAndPrediction_H_
#define ConvolutionAndPrediction_H_

void getDiffusionSignal(double **q, double **v, int qLen, int vLen,double**out); //IMPORTANT: Output
// Matrix is dynamically allocated, must free up when done

double* estimateWeightVector(double ** H, int qLen, int vLen, double * e, int eLen);
//Also dynamically allocated

void getReconstructionMatrix(double ** G, double ** H, int k, int n, int p, int m, double ** out);

double calculateNormalizationConstant(double ** A, double * e, int n, int m);

void computeODF(double ** A, double * e, int n, int m, double * out);

#endif
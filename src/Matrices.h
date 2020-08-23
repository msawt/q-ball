
double mulPoints(double* rowM, double* colM, int length);

void matMulPoints(double p1[3], double p2[3], double out[3][3]);

void div33ByScalar(double m1[3][3], double div, double out[3][3]);

void mulThreeByThree(double m1[3][3], double m2[3][3], double out[3][3]);

void pseudoInverse(double ** mat, int rowDim, int colDim, double ** out);

void ADJ(double ** M,double ** adj, int N);

void getCfactor(double ** M, double ** t, int n);

double det(double ** A, int n);

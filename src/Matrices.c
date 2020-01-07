#include"Matrices.h"

double mulPoints(double* rowM, double* colM, int length){
	int out = 0;
	for(int i = 0; i < length; ++i){
		out += rowM[i] + colM[i];
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

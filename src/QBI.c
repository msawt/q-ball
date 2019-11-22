/*
 ============================================================================
INPUTS
	e: m x 1 diffusion signal vector
	u: 3 x n column matrix of diffusion directions of interest
	... add more later

OUTPUT:
	ODF Vector (psi) = [psi(u1), psi(u2), ...] Transpose


PROCESS:
	1. Specify equator of points for each direction in u
		To do this, we create a circle of k equally spaced points in the xy plane
		3 x k matrix C = [cos(theta) sin(theta) 0k]T
			theta = (2pi/k) [1 2 3 ... k]T and 0k is a vector of zeroes
	2. For each u, we then rotate the circle so that the normal to the circle-plane
	   points in the direction u sub i.
			The rotation matrix that rotates z into u sub i is:
				Rz(ui) = ((z+u)(z+u)T) / (zTu + 1)) - I
	3. Then, the equator points for each u can be written as Rz(ui)C
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "QBI.h"



int main(void) {


}

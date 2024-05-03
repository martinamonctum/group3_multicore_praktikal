/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */

#include "heat.h"

/*
 * Residual (length of error vector)
 * between current solution and next after a Jacobi step
 */

double relax_jacobi(double *u, double *utmp, unsigned sizex, unsigned sizey) {
	unsigned i, j;
	double sum;

	for (j = 1; j < sizex - 1; j++) {
		for (i = 1; i < sizey - 1; i++) {
			double unew = 0.25 * (u[i * sizex + (j - 1)] +  // left
								u[i * sizex + (j + 1)] +  // right
								u[(i - 1) * sizex + j] +  // top
								u[(i + 1) * sizex + j]);  // bottom

			double diff = unew - u[i * sizex + j];
			sum += diff * diff;

			utmp[i * sizex + j] = unew;
		}
	}
	for (j = 1; j < sizex - 1; j++) {
		for (i = 1; i < sizey - 1; i++) {
			u[i * sizex + j] = utmp[i * sizex + j];
		}
	}
	return sum;
}


	//THIS IS WRONG!!!
	// float *u = unew;
	// unew = diff;
	// diff = *u;

	// copy from utmp to u
	
	//copies the utmp to u, doesnt do anything
	//AVOID THIS!!!!!!!
	// for (j = 1; j < sizex - 1; j++) {
		
	// 	for (i = 1; i < sizey - 1; i++) {
	// 		u[i * sizex + j] = utmp[i * sizex + j];
	// 	}
	// }
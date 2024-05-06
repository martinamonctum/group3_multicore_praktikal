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
double residual_jacobi(double *u, unsigned sizex, unsigned sizey) {
	unsigned i, j;
	double unew, diff, sum = 0.0;

	for (j = 1; j < sizex - 1; j++) {
		for (i = 1; i < sizey - 1; i++) {
			unew = 0.25 * (u[i * sizex + (j - 1)] +  // left
						u[i * sizex + (j + 1)] +  // right
						u[(i - 1) * sizex + j] +  // top
						u[(i + 1) * sizex + j]); // bottom

			diff = unew - u[i * sizex + j];
			sum += diff * diff;
		}
	}

	return sum;
}

/*
 * One Jacobi iteration step
 */
void relax_jacobi(double *u, double *utmp,unsigned sizex, unsigned sizey) {

	int i, j;

    // this only works for even resólution, for uneven resolution the loops have to be split, due to the borders.
	for (i = 1; i < sizey - 1; i+=4) {
	    for (j = 2; j < sizex - 1; j+=2) {
			u[i * sizex + j] = 0.25 * (u[i * sizex + (j - 1)] +  // left
						u[i * sizex + (j + 1)] +  // right
						u[(i - 1) * sizex + j] +  // top
						u[(i + 1) * sizex + j]); // bottom
            u[(i+1) * sizex + (j-1)] = 0.25 * (u[(i+1) * sizex + (j-2)] +  // left
						u[(i+1) * sizex + j] +  // right
						u[i * sizex + (j-1)] +  // top
						u[(i+2) * sizex + (j-1)]); // bottom
		}
        for (j = sizex-2; j > 0; j-=2) {
			u[(i+2) * sizex + j] = 0.25 * (u[(i+2) * sizex + (j - 1)] +  // left
						u[(i+2) * sizex + (j + 1)] +  // right
						u[(i+1) * sizex + j] +  // top
						u[(i+3) * sizex + j]); // bottom
            u[(i+3) * sizex + (j-1)] = 0.25 * (u[(i+3) * sizex + (j-2)] +  // left
						u[(i+3) * sizex + j] +  // right
						u[(i+2) * sizex + (j-1)] +  // top
						u[(i+4) * sizex + (j-1)]); // bottom
		}
	}

    for (i = 1; i < sizey - 1; i+=4) {
	    for (j = 1; j < sizex - 1; j+=2) {
            u[i * sizex + j] = 0.25 * (u[i * sizex + (j - 1)] +  // left
						u[i * sizex + (j + 1)] +  // right
						u[(i - 1) * sizex + j] +  // top
						u[(i + 1) * sizex + j]); // bottom
			u[(i+1) * sizex + (j+1)] = 0.25 * (u[(i+1) * sizex + j] +  // left
						u[(i+1) * sizex + (j + 2)] +  // right
						u[i * sizex + (j+1)] +  // top
						u[(i+2) * sizex + (j+1)]); // bottom
            
		}
        for (j = sizex-3; j > 0; j-=2) {
            u[(i+2) * sizex + j] = 0.25 * (u[(i+2) * sizex + (j - 1)] +  // left
						u[(i+2) * sizex + (j + 1)] +  // right
						u[(i+1) * sizex + j] +  // top
						u[(i+3) * sizex + j]); // bottom
			u[(i+3) * sizex + (j+1)] = 0.25 * (u[(i+3) * sizex + j] +  // left
						u[(i+3) * sizex + (j + 2)] +  // right
						u[(i+2) * sizex + (j+1)] +  // top
						u[(i+4) * sizex + (j+1)]); // bottom
            
		}
	}

}
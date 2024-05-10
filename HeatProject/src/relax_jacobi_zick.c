/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */

#include "heat.h"

/*
 * One Jacobi iteration step
 */
double relax_jacobi(double **u1, double **utmp1,
         unsigned sizex, unsigned sizey) {

	int i, j;
    double *help,*u, *utmp,factor=0.5;

    utmp=*utmp1;
    u=*u1;
    double unew, diff, sum=0.0;

    // this only works for even resólution, for uneven resolution the loops have to be split, due to the borders.
	for (i = 1; i < sizey - 1; i+=2) {
        int ii=i*sizex;
  	    int iim1=(i-1)*sizex;
  	    int iip1=(i+1)*sizex;
        int iip2=(i+2)*sizex;
		#pragma ivdep
	    for (j = 2; j < sizex - 1; j+=2) {
			unew = 0.25 * (u[ ii+(j-1) ]+
        		            u[ ii+(j+1) ]+
        		            u[ iim1+j ]+
        		            u[ iip1+j ]);
		    diff = unew - u[ii + j];
		    utmp[ii+j] = unew;
		    sum += diff * diff;

            unew2 = 0.25 * (u[ iip1+(j-2) ]+
        		            u[ iip1+(j) ]+
        		            u[ ii+(j-1) ]+
        		            u[ iip2+(j-1) ]);
		    diff2 = unew2 - u[iip1 + (j-1)];
		    utmp[iip1+(j-1)] = unew2;
		    sum += diff2 * diff2;
		}
	}

    for (i = 1; i < sizey - 1; i+=2) {
        int ii=i*sizex;
  	    int iim1=(i-1)*sizex;
  	    int iip1=(i+1)*sizex;
        int iip2=(i+2)*sizex;
		#pragma ivdep
	    for (j = 1; j < sizex - 1; j+=2) {
            unew = 0.25 * (u[ ii+(j-1) ]+
        		            u[ ii+(j+1) ]+
        		            u[ iim1+j ]+
        		            u[ iip1+j ]);
		    diff = unew - u[ii + j];
		    utmp[ii+j] = unew;
		    sum += diff * diff;

			unew2 = 0.25 * (u[ iip1+(j) ]+
        		            u[ iip1+(j+2) ]+
        		            u[ ii+(j+1) ]+
        		            u[ iip2+(j+1) ]);
		    diff2 = unew2 - u[ iip1+(j+1) ];
		    utmp[iip1+(j+1)] = unew2;
		    sum += diff2 * diff2;
            
		}
	}
	*u1=utmp;
  	*utmp1=u;
	return(sum);
}
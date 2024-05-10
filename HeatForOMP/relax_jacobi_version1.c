/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */

#include "heat.h"
#include <omp.h>


double relax_jacobi( double **u1, double **utmp1,
         unsigned sizex, unsigned sizey )
{

  omp_set_num_threads(4);
  //set number of threads for a parallel region
  
  int i, j, ii, iim1, iip1;
  double *help,*u, *utmp,factor=0.5, sum=0.0;
  double unew, diff;
  utmp=*utmp1;
  u=*u1;


  //#pragma omp parallel{}
  #pragma omp parallel for schedule(static) reduction(+: sum) private(i, j, ii, iim1, iip1, unew, diff)

  for( int i=1; i<sizey-1; i++ ) {
    ii=i*sizex;
    iim1=(i-1)*sizex;
    iip1=(i+1)*sizex;
    
    
    for(int j=1; j<sizex-1; j++ ){
      
      unew = 0.25 * (u[ ii+(j-1) ]+
                        u[ ii+(j+1) ]+
                        u[ iim1+j ]+
                        u[ iip1+j ]);
        diff = unew - u[ii + j];
        utmp[ii+j] = unew;
        sum += diff * diff;

    }
  }
  

  *u1=utmp;
  *utmp1=u;
  return(sum);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////





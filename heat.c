/*
 * heat.h
 *
 * Iterative solver for heat distribution
 */

#include "heat.h"

#include "papi.h"
#include <stdio.h>
#include <stdlib.h>

#include "input.h"
#include "timing.h"

#define NUM_FLOPS 10000
#define NUM_EVENTS 1

void usage(char *s) {
	fprintf(stderr, "Usage: %s <input file> [result file]\n\n", s);
}

int main(int argc, char *argv[]) {
	unsigned iter;
	FILE *infile, *resfile;
	char *resfilename;

	// algorithmic parameters
	algoparam_t param;
	int np,i;

	double runtime, flop;
	double residual;
	double time[1000];
	double floprate[1000];
	double papi_floprate[1000];
	int resolution[1000];
	int experiment=0;

	// check arguments
	if (argc < 2) {
		usage(argv[0]);
		return 1;
	}

	// check input file
	if (!(infile = fopen(argv[1], "r"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);

		usage(argv[0]);
		return 1;
	}

	// check result file
	resfilename = (argc >= 3) ? argv[2] : "heat.ppm";

	if (!(resfile = fopen(resfilename, "w"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for writing.\n\n", resfilename);

		usage(argv[0]);
		return 1;
	}

	// check input
	if (!read_input(infile, &param)) {
		fprintf(stderr, "\nError: Error parsing input file.\n\n");

		usage(argv[0]);
		return 1;
	}

	print_params(&param);

	// set the visualization resolution
	param.visres = 1024;

	param.u = 0;
	param.uhelp = 0;
	param.uvis = 0;

	param.act_res = param.initial_res;
	char c = 0;
	// loop over different resolutions
	while (1) {

		// free allocated memory of previous experiment
		if (param.u != 0)
			finalize(&param);

		if (!initialize(&param)) {
			fprintf(stderr, "Error in Jacobi initialization.\n\n");

			usage(argv[0]);
		}

		fprintf(stderr, "Resolution: %5u\r", param.act_res);

		// full size (param.act_res are only the inner points)
		np = param.act_res + 2;

		// starting time
		runtime = wtime();
		residual = 999999999;

		int retval;
		retval = PAPI_hl_region_begin("computation");
		if (retval!= PAPI_OK){
			printf("PAPI IS NOT OK!!!!");
			return 0;
		}
		

		//float real_time, proc_time, mflops;
  		// long long flpops;
		// int flops_retval;
		
		//flops_retval = PAPI_flops_rate(PAPI_FP_OPS, &real_time, &proc_time, &flpops, &mflops);
		// flops_retval = mflops;
		// if(flops_retval < PAPI_OK){
		// 	printf("PAPI IS NOT OK!!!!");
		// 	return 0;
		// }
		iter = 0;
		while (1) {
			float mflops = PAPI_DP_OPS;
			switch (param.algorithm) {

			case 0: // JACOBI

				relax_jacobi(param.u, param.uhelp, np, np);
				residual = residual_jacobi(param.u, np, np);
				break;

			case 1: // GAUSS

				relax_gauss(param.u, np, np);
				residual = residual_gauss(param.u, param.uhelp, np, np);
				break;
			}

			iter++;

			// solution good enough ?
			if (residual < 0.000005)
				break;

			// max. iteration reached ? (no limit with maxiter=0)
			if (param.maxiter > 0 && iter >= param.maxiter)
				break;

			// if (iter % 100 == 0)
			// 	fprintf(stderr, "residual %f, %d iterations\n", residual, iter);
		}
		//flops_retval = PAPI_flops_rate(PAPI_FP_OPS, &real_time, &proc_time, &flpops, &mflops);
		// flops_retval = mflops;
		// if(flops_retval < PAPI_OK){
		// 	printf("PAPI IS NOT OK!!!!");
		// 	return 0;
		// }
		float mflops;
		
		
		

		// Flop count after <i> iterations
		flop = iter * 11.0 * param.act_res * param.act_res;
		// stopping time
		runtime = wtime() - runtime;
		mflops = PAPI_DP_OPS/(PAPI_TOT_CYC*1000000); //(PAPI_TOT_CYC*1000000);
		
		papi_floprate[experiment] = mflops;
		retval = PAPI_hl_region_end("computation");
		if (retval!= PAPI_OK){
			printf("PAPI IS NOT OK!!!!");
			return 0;
		}

		fprintf(stderr, "Resolution: %5u, ", param.act_res);
		fprintf(stderr, "Time: %04.3f ", runtime);
		fprintf(stderr, "(%3.3f GFlop => %6.2f MFlop/s, ", flop / 1000000000.0, flop / runtime / 1000000);
		fprintf(stderr, "residual %f, %d iterations)\n", residual, iter);

		// for plot...
		time[experiment]=runtime;
		floprate[experiment]=flop / runtime / 1000000;
		resolution[experiment]=param.act_res;
		experiment++;
		c++;
		if (param.act_res + param.res_step_size > param.max_res)
			break;
		if (c==5)
			break;
		param.act_res += param.res_step_size;
		
	}

    FILE *flopDat = fopen("FlopData", "w");

	for (i=0;i<experiment; i++){
		fprintf(flopDat, "%5d " ,resolution[i]);
        fprintf(flopDat, "%5.3f ", floprate[i]);
		fprintf(flopDat, "%5.3f\n", papi_floprate[i]);
		printf("%5d; %5.3f; %5.3f\n", resolution[i], time[i], floprate[i]);
	}
	fclose(flopDat);
	
	coarsen(param.u, np, np, param.uvis, param.visres + 2, param.visres + 2);
	write_image(resfile, param.uvis, param.visres + 2, param.visres + 2);

	finalize(&param);

      return 0;
}

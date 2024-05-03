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
	// // papi init
	// int retval,status=0, EventSet=PAPI_NULL;
	// int EventCode[4] = {PAPI_L2_TCA,PAPI_L2_TCM,PAPI_L3_TCA,PAPI_L3_TCM};
	// long_long values[4];
	// retval = PAPI_library_init(PAPI_VER_CURRENT);
	// if (retval != PAPI_VER_CURRENT) {
	// fprintf(stderr, "PAPI library init error!\n");
	// exit(1);
	// }
	// /* Create the Event Set */
	// if (PAPI_create_eventset(&EventSet) != PAPI_OK){
	// 	fprintf(stderr, "PAPI create event error\n");
	// 	exit(1);
	// }
	// /* Add Total Instructions Executed to our Event Set */
	// if (PAPI_add_events(EventSet, EventCode,4) != PAPI_OK){
	// 	fprintf(stderr, "PAPI add event error!\n");
	// 	exit(1);
	// }
	// set the visualization resolution
	param.visres = 1024;

	param.u = 0;
	param.uhelp = 0;
	param.uvis = 0;

	param.act_res = param.initial_res;
	char c = 0;
	// if (PAPI_start(EventSet) != PAPI_OK){
	// 		fprintf(stderr, "PAPI start event error!\n");
	// 		exit(1);
	// 	}
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
		/* Start counting events in the Event Set */
		
		double ftime;
		ftime = PAPI_get_real_usec();
		
		iter = 0;
		while (1) {

			switch (param.algorithm) {

			case 0: // JACOBI

				residual = relax_jacobi(param.u, param.uhelp, np, np);
				//residual = residual_jacobi(param.u, np, np);
				printf("REDISUAL:  %f", residual);
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
		retval = PAPI_hl_region_end("computation");
		if (retval!= PAPI_OK){
			printf("PAPI IS NOT OK!!!!");
			return 0;
		}
		/* Read the counting events in the Event Set */
		// if (PAPI_read(EventSet, values) != PAPI_OK){
		// 	fprintf(stderr, "PAPI add event error!\n");
		// 	exit(1);
		// }
		// FILE *cacheMisses = fopen("cacheMisses", "a");
		// for (i=0;i<4; i++){
		// 	fprintf(cacheMisses, "%5d " ,values[i]);
		// }
		// fprintf(cacheMisses,"\n");
		// fprintf(cacheMisses, "L2 miss rate%5d\n", values[1]/values[0]);
		// fprintf(cacheMisses, "L3 miss rate%5d\n", values[3]/values[2]);
		
		ftime = PAPI_get_real_usec() - ftime;
		// papi_floprate[experiment] = values[0]/ftime;
		
		// Flop count after <i> iterations
		flop = iter * 11.0 * param.act_res * param.act_res;
		// stopping time
		runtime = wtime() - runtime;

		fprintf(stderr, "Resolution: %5u, ", param.act_res);
		fprintf(stderr, "Time in microseconds: %04.3f ", ftime);
		fprintf(stderr, "(%3.3f GFlop => %6.2f MFlop/s, ", flop / 1000000.0, flop / runtime / 1000000);
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
		
		// /* Reset the counting events in the Event Set */
		// if (PAPI_reset(EventSet) != PAPI_OK){
		// 	fprintf(stderr, "PAPI reset event error!\n");
		// 	exit(1);
		// }
	}
	
    FILE *flopDat = fopen("FlopData", "w");

	for (i=0;i<experiment; i++){
		fprintf(flopDat, "%5d " ,resolution[i]);
        fprintf(flopDat, "%5.3f ", floprate[i]);
		fprintf(flopDat, "%5.3f\n", papi_floprate[i]);
		printf("%5d; %5.3f; %5.3f; %5.3f\n", resolution[i], time[i], floprate[i], papi_floprate[i]);
	}
	fclose(flopDat);
	
	coarsen(param.u, np, np, param.uvis, param.visres + 2, param.visres + 2);
	write_image(resfile, param.uvis, param.visres + 2, param.visres + 2);

	finalize(&param);
	return 0;
}

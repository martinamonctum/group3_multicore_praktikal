#include <stdio.h>
#include <stdlib.h>

#include "input.h"
#include "heat.h"
#include "timing.h"

double* time;

void usage(char *s) {
	fprintf(stderr, "Usage: %s <input file> [result file]\n\n", s);
}

int main(int argc, char *argv[]) {

	int err;

	int rank;
	int coords[2];
	int size;
	MPI_Status status;

	MPI_Comm comm_cart;
	int dims[2];
	int periods[2];
	int xdim;
	int ydim;
	int xshift;
	int yshift;

	int left_rank;
	int right_rank;
	int upper_rank;
	int lower_rank;

	int left_pos[2];
	int right_pos[2];
	int upper_pos[2];
	int lower_pos[2];

	double* left_border_send;
	double* right_border_send;
	double* lower_border_send;
	double* upper_border_send;

	double* left_border_rec;
	double* right_border_rec;
	double* lower_border_rec;
	double* upper_border_rec;

	int i, j, k, m, n, ret;
	FILE *infile, *resfile;
	char *resfilename;
	int np, iter, chkflag;
	
	// algorithmic parameters
	algoparam_t param;

	// timing

	double residual;
	double global_residual;

	// set the visualization resolution
	param.visres = 100;

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
	time = (double *) calloc(sizeof(double), (int) (param.max_res - param.initial_res + param.res_step_size) / param.res_step_size);

	
	int exp_number = 0;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	dims[0] = 2;
	dims[1] = 2;

	periods[0] = 0;
	periods[1] = 0;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_cart);

	MPI_Cart_coords(comm_cart, rank, 2, coords);

	left_pos[0] = coords[0]; 
	left_pos[1] = coords[1]-1;
	right_pos[0] = coords[0];
	right_pos[1] = coords[1]+1;
	upper_pos[0] = coords[0]-1;
	upper_pos[1] = coords[1];
	lower_pos[0] = coords[0]+1;
	lower_pos[1] = coords[1];

	// Returns MPI_PROC_NULL for wrong coords, which is fine.
	if (left_pos[1] >= 0) MPI_Cart_rank(comm_cart, left_pos, &left_rank); else left_rank = MPI_PROC_NULL;
	if (right_pos[1] < dims[1]) MPI_Cart_rank(comm_cart, right_pos, &right_rank); else right_rank = MPI_PROC_NULL;
	if (upper_pos[0] >= 0) MPI_Cart_rank(comm_cart, upper_pos, &upper_rank); else upper_rank = MPI_PROC_NULL;
	if (lower_pos[0] < dims[0]) MPI_Cart_rank(comm_cart, lower_pos, &lower_rank); else lower_rank = MPI_PROC_NULL;

	printf("I am %d! And live at (%d,%d)\n", rank, coords[0], coords[1]);
	

	for (param.act_res = param.initial_res; param.act_res <= param.max_res; param.act_res = param.act_res + param.res_step_size) {
		if (param.act_res!=param.initial_res) finalize(&param);

		np = param.act_res + 2;

		if (coords[1] < dims[1] - 1){
			xdim = param.act_res / dims[1];
		} else {
			xdim = param.act_res - param.act_res / dims[1];
		}
		xshift = param.act_res / dims[1];
		xdim += 2;

		if (coords[0] < dims[0] - 1){
			ydim = param.act_res /dims[0];
		} else {
			ydim = param.act_res - param.act_res / dims[0];
		}
		ydim += 2;
		yshift = param.act_res / dims[0];

		if (!initialize(&param,xdim,ydim,xshift,yshift,coords,dims)) {
			fprintf(stderr, "Error in Jacobi initialization.\n\n");

			usage(argv[0]);
		}

		for (i = 0; i < ydim; i++) {
			for (j = 0; j < xdim; j++) {
				param.uhelp[i * xdim + j] = param.u[i * xdim + j];
			}
		}

		lower_border_send = (double*)malloc( sizeof(double)* (xdim-2) );
		lower_border_rec = (double*)malloc( sizeof(double)* (xdim-2) );
		upper_border_send = (double*)malloc( sizeof(double)* (xdim-2) );
		upper_border_rec = (double*)malloc( sizeof(double)* (xdim-2) );
		left_border_send = (double*)malloc( sizeof(double)* (ydim-2) );
		left_border_rec = (double*)malloc( sizeof(double)* (ydim-2) );
		right_border_send = (double*)malloc( sizeof(double)* (ydim-2) );
		right_border_rec = (double*)malloc( sizeof(double)* (ydim-2) );

		// starting time
		time[exp_number] = MPI_Wtime();
		residual = 999999999;

		for (iter = 0; iter < param.maxiter; iter++) {
			global_residual = 0;
			residual = relax_jacobi(&(param.u), &(param.uhelp), xdim, ydim);
			MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
			if(rank==0){
				printf("residual iter %d: %f\n", iter, global_residual);
			}
			if (global_residual<0.00000005)break;
			for(m = 1; m < xdim-1; m++){
				// send not border but newly computed data.
			    lower_border_send[m] = param.u[(ydim-2)*xdim+m];
				upper_border_send[m] = param.u[xdim+m];
			}
			for(n = 1; n < ydim-1; n++){
				// send not border but newly computed data.
				left_border_send[n] = param.u[n*xdim+1];
				right_border_send[n] = param.u[n*xdim+(xdim-2)];
			}

			MPI_Sendrecv(lower_border_send, xdim-2, MPI_DOUBLE, lower_rank, 1,
						lower_border_rec, xdim-2, MPI_DOUBLE, lower_rank, 1,
						comm_cart, MPI_STATUS_IGNORE);

			MPI_Sendrecv(upper_border_send, xdim-2, MPI_DOUBLE, upper_rank, 1,
						upper_border_rec, xdim-2, MPI_DOUBLE, upper_rank, 1,
						comm_cart, MPI_STATUS_IGNORE);

			MPI_Sendrecv(left_border_send, ydim-2, MPI_DOUBLE, left_rank, 1,
						left_border_rec, ydim-2, MPI_DOUBLE, left_rank, 1,
						comm_cart, MPI_STATUS_IGNORE);

			MPI_Sendrecv(right_border_send, ydim-2, MPI_DOUBLE, right_rank, 1,
						right_border_rec, ydim-2, MPI_DOUBLE, right_rank, 1,
						comm_cart, MPI_STATUS_IGNORE);

			for(m = 1; m < xdim-1; m++){
				if (lower_rank != MPI_PROC_NULL) param.u[(ydim-1)*xdim+m] = lower_border_rec[m];
				if (upper_rank != MPI_PROC_NULL) param.u[m] = upper_border_rec[m];
			}
			for(n = 1; n < ydim-1; n++){
				if (left_rank != MPI_PROC_NULL) param.u[n*xdim] = left_border_rec[n];
				if (right_rank != MPI_PROC_NULL) param.u[n*xdim+(xdim-1)] = right_border_rec[n];
			}

			MPI_Barrier(comm_cart);

		}

		time[exp_number] = MPI_Wtime() - time[exp_number];

		if (rank == 0){
			printf("\n\nResolution: %u\n", param.act_res);
			printf("===================\n");
			printf("Execution time: %f\n", time[exp_number]);
			printf("Residual: %f\n\n", global_residual);

			printf("megaflops:  %.1lf\n", (double) param.maxiter * (np - 2) * (np - 2) * 7 / time[exp_number] / 1000000);
			printf("  flop instructions (M):  %.3lf\n", (double) param.maxiter * (np - 2) * (np - 2) * 7 / 1000000);
		}

		free(lower_border_send);
		free(lower_border_rec);
		free(upper_border_send);
		free(upper_border_rec);
		free(left_border_send);
		free(left_border_rec);
		free(right_border_send);
		free(right_border_rec);

		exp_number++;
	}

	err = MPI_Finalize();

    param.uvis  = (double*)calloc( sizeof(double),
				      (param.visres+2) *
				      (param.visres+2) );

	param.act_res = param.act_res - param.res_step_size;

	//coarsen(param.u, param.act_res + 2, param.act_res + 2, param.uvis, param.visres + 2, param.visres + 2);

	//write_image(resfile, param.uvis, param.visres + 2, param.visres + 2);


	return 0;
}

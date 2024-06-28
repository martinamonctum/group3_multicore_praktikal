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
	MPI_Request req;

	MPI_Request reqs[9];

	MPI_Comm comm_cart;
	int dims[2];
	int periods[2];
	int xdim;
	int ydim;
	int xshift;
	int yshift;

	int xdim_vis;
	int ydim_vis;
	int xshift_vis;
	int yshift_vis;

	int vis_partner_pos[2];
	int vis_partner_rank;

	int local_xpos_vis;
	int local_ypos_vis;

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

	int csv_flag = 0;

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

	time = (double *) calloc(sizeof(double), (int) (param.max_res - param.initial_res + param.res_step_size) / param.res_step_size);

	
	int exp_number = 0;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc >= 5){
		dims[0] = atoi(argv[3]);
		dims[1] = atoi(argv[4]);
	} else {
		fprintf(stderr, "\nError: This is the MPI version! You need to give the Cartesian Grid dimensions!\n\n");
		return 1;
	}

	if (argc >=6 && argv[5][0] == 'c' && argv[5][1] == 's' && argv[5][2] == 'v') csv_flag = 1;

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
	
	if (rank == 0 && !csv_flag) print_params(&param);

	for (param.act_res = param.initial_res; param.act_res <= param.max_res; param.act_res = param.act_res + param.res_step_size) {
		if (param.act_res!=param.initial_res) finalize(&param);

		np = param.act_res + 2;

		if (coords[1] < dims[1] - 1){
			xdim = param.act_res / dims[1];
			xdim_vis = param.visres / dims[1];
		} else {
			xdim = param.act_res - (dims[1]-1)*(param.act_res / dims[1]);
			xdim_vis = param.visres - (dims[1]-1)*(param.visres / dims[1]);
		}
		xshift = param.act_res / dims[1];
		xdim += 2;
		xdim_vis += 2;
		xshift_vis = param.visres / dims[1];

		if (coords[0] < dims[0] - 1){
			ydim = param.act_res /dims[0];
			ydim_vis = param.visres / dims[0];
		} else {
			ydim = param.act_res - (dims[0]-1)*(param.act_res / dims[0]);
			ydim_vis = param.visres - (dims[0]-1)*(param.visres / dims[0]);
		}
		ydim += 2;
		yshift = param.act_res / dims[0];
		ydim_vis += 2;
		yshift_vis = param.visres / dims[0];

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
			
			for(m = 1; m < xdim-1; m++){
				// send not border but newly computed data.
			    lower_border_send[m-1] = param.u[(ydim-2)*xdim+m];
				upper_border_send[m-1] = param.u[xdim+m];
			}
			for(n = 1; n < ydim-1; n++){
				// send not border but newly computed data.
				left_border_send[n-1] = param.u[n*xdim+1];
				right_border_send[n-1] = param.u[n*xdim+(xdim-2)];
			}

			MPI_Irecv(lower_border_rec, xdim-2, MPI_DOUBLE, lower_rank, 1,
					comm_cart, &reqs[0]);

			MPI_Isend(lower_border_send, xdim-2, MPI_DOUBLE, lower_rank, 3,
					comm_cart, &reqs[1]);

			MPI_Irecv(upper_border_rec, xdim-2, MPI_DOUBLE, upper_rank, 3,
					comm_cart, &reqs[2]);

			MPI_Isend(upper_border_send, xdim-2, MPI_DOUBLE, upper_rank, 1,
					comm_cart, &reqs[3]);

			MPI_Irecv(left_border_rec, ydim-2, MPI_DOUBLE, left_rank, 2,
					comm_cart, &reqs[4]);

			MPI_Isend(left_border_send, ydim-2, MPI_DOUBLE, left_rank, 0,
					comm_cart, &reqs[5]);

			MPI_Irecv(right_border_rec, ydim-2, MPI_DOUBLE, right_rank, 0,
					comm_cart, &reqs[6]);

			MPI_Isend(right_border_send, ydim-2, MPI_DOUBLE, right_rank, 2,
					comm_cart, &reqs[7]);
			

			MPI_Iallreduce(&residual, &global_residual, 1, MPI_DOUBLE, MPI_SUM, comm_cart, &reqs[8]);

			MPI_Waitall(9, reqs, MPI_STATUSES_IGNORE);

			if (global_residual<0.00000005)break;

			for(m = 1; m < xdim-1; m++){
				if (lower_rank != MPI_PROC_NULL) param.u[(ydim-1)*xdim+m] = lower_border_rec[m-1];
				if (upper_rank != MPI_PROC_NULL) param.u[m] = upper_border_rec[m-1];
			}
			for(n = 1; n < ydim-1; n++){
				if (left_rank != MPI_PROC_NULL) param.u[n*xdim] = left_border_rec[n-1];
				if (right_rank != MPI_PROC_NULL) param.u[n*xdim+(xdim-1)] = right_border_rec[n-1];
			}

		}

		time[exp_number] = MPI_Wtime() - time[exp_number];

		if (rank == 0){
			if (!csv_flag){
				printf("\n\nResolution: %u\n", param.act_res);
				printf("===================\n");
				printf("Execution time: %f\n", time[exp_number]);
				printf("Residual: %f\n\n", global_residual);

				printf("megaflops:  %.1lf\n", (double) param.maxiter * (np - 2) * (np - 2) * 7 / time[exp_number] / 1000000);
				printf("  flop instructions (M):  %.3lf\n", (double) param.maxiter * (np - 2) * (np - 2) * 7 / 1000000);
			} else {
				printf("%d,%d,%f\n",size,param.act_res,time[exp_number]);
			}
		} 

		exp_number++;

		free(lower_border_send);
		free(upper_border_send);
		free(left_border_send);
		free(right_border_send);
		free(lower_border_rec);
		free(upper_border_rec);
		free(left_border_rec);
		free(right_border_rec);
		lower_border_send=0;
		upper_border_send=0;
		left_border_send=0;
		right_border_send=0;
		lower_border_rec=0;
		upper_border_rec=0;
		left_border_rec=0;
		right_border_rec=0;
	}

    param.uvis  = (double*)calloc( sizeof(double),
				      ydim_vis *
				      xdim_vis );

	coarsen(param.u, xdim, ydim, param.uvis, xdim_vis, ydim_vis);

	// This is annoying, actually only rank 0 would need this.
	// Hope is: exchanging calloc for malloc makes it such that only rank 0 actually allocates.
	double *uvis_global = (double*)malloc( sizeof(double)*
				      (param.visres+2) *
				      (param.visres+2) );

	MPI_Request *img_reqs = (MPI_Request*) malloc(sizeof(MPI_Request)*
						(param.visres+2) *
						(param.visres+2));
	if (rank == 0){

		for (i = 0; i < param.visres+2; i++){
			for(j = 0; j < param.visres+2; j++){
				vis_partner_pos[0] = i / yshift_vis;
				vis_partner_pos[1] = j / xshift_vis;
				if (vis_partner_pos[0] > dims[0]-1) vis_partner_pos[0] = dims[0]-1;
				if (vis_partner_pos[1] > dims[1]-1) vis_partner_pos[1] = dims[1]-1;
				MPI_Cart_rank(comm_cart, vis_partner_pos, &vis_partner_rank);
				MPI_Irecv(&uvis_global[i*(param.visres+2)+j], 1, MPI_DOUBLE, vis_partner_rank, (i)*(param.visres+2)+(j), comm_cart, &img_reqs[i*(param.visres+2)+j]);
			}
		}
	}

	for (i = 0; i < param.visres+2; i++){
		for(j = 0; j < param.visres+2; j++){
			vis_partner_pos[0] = i / yshift_vis;
			vis_partner_pos[1] = j / xshift_vis;
			if (vis_partner_pos[0] > dims[0]-1) vis_partner_pos[0] = dims[0]-1;
			if (vis_partner_pos[1] > dims[1]-1) vis_partner_pos[1] = dims[1]-1;
			MPI_Cart_rank(comm_cart, vis_partner_pos, &vis_partner_rank);
			if (rank == vis_partner_rank){
				local_ypos_vis = i - vis_partner_pos[0] * yshift_vis;
				local_xpos_vis = j - vis_partner_pos[1] * xshift_vis;
				MPI_Isend(&param.uvis[local_ypos_vis*xdim_vis+local_xpos_vis], 1, MPI_DOUBLE, 0, i*(param.visres+2)+j, comm_cart, &req);
			}
		}
	}

	if (rank==0) MPI_Waitall((param.visres+2)*(param.visres+2), img_reqs, MPI_STATUSES_IGNORE);

	if (rank == 0) write_image(resfile, uvis_global, param.visres+2, param.visres+2);

	free(uvis_global);
	uvis_global=0;

	free(param.uvis);
	param.uvis=0;

	finalize(&param);

	err = MPI_Finalize();

	return 0;
}

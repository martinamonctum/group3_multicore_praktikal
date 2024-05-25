#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    int num_procs;
    int rank;
    int n_iter = 25;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int s_buffer[n_iter];
    for(int i=0;i<n_iter;i++) *(s_buffer+i) = (rank+1)*i;
    int r_buffer[n_iter];

    double time = 0.;
    if (rank == 0) {
        int dest = 1;
        int source = 1;

        MPI_Barrier(MPI_COMM_WORLD);

        time -= MPI_Wtime();
        for (int i = 0; i < n_iter; i++) {
            MPI_Send(s_buffer+i,
                                 1, MPI_INT, dest,
                                 0, MPI_COMM_WORLD);

            MPI_Recv(r_buffer+i,
                                 1, MPI_INT, source,
                                 0, MPI_COMM_WORLD, &status);

        } 
        time += MPI_Wtime();
    } else if (rank == 1) {
        int dest = 0;
        int source = 0;

        MPI_Barrier(MPI_COMM_WORLD);

        time -= MPI_Wtime();
        for (int i = 0; i < n_iter; i++) {
            MPI_Recv(r_buffer+i,
                                 1, MPI_INT, source,
                                 0, MPI_COMM_WORLD, &status);

            MPI_Send(s_buffer+i,
                                 1, MPI_INT, dest,
                                 0, MPI_COMM_WORLD);

        }
        time += MPI_Wtime();
    }
    time /= n_iter;

    MPI_Finalize();

    printf("I am %d, took: %f\n", rank, time);
    printf("Result: ");
    for (int i = 0; i < n_iter; i++){
        printf("%d ", r_buffer[i]);
    }
    printf("\n");

}
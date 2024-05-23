#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
   const int PING_PONG_LIMIT = 10;
   
   // Initialize the MPI environment
   MPI_Init(NULL, NULL);
   // Find out rank, size
   int process_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
   int process_num;
   MPI_Comm_size(MPI_COMM_WORLD, &process_num);

   // We are assuming 2 processes for this task
   if (process_num != 2) {
      fprintf(stderr, "World size must be two for %s\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   int ping_pong_count = 0;
   int partner_rank = (process_rank + 1) % 2;
   MPI_Barrier(MPI_COMM_WORLD);
   double start,time;
   while (ping_pong_count < PING_PONG_LIMIT) {
      if (process_rank == ping_pong_count % 2) {
         start = MPI_Wtime();
         // Increment the ping pong count before you send it
         ping_pong_count++;

         MPI_Send(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
         time += MPI_Wtime() - start;
         printf("process %d sent and incremented ping_pong_count %d\n",
               process_rank, ping_pong_count);
      } else {
         start = MPI_Wtime();
         MPI_Recv(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD,
                  MPI_STATUS_IGNORE);
         time += MPI_Wtime() - start;
         printf("process %d received ping_pong_count %d\n",
               process_rank, ping_pong_count);
      }
   }
   MPI_Finalize();
   printf("time of sending the message back and forth is %d ", time / 10);
   printf("my rank is:%d\n",process_rank);
}
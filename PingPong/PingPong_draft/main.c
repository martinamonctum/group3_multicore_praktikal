#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
   const int PING_PONG_LIMIT = 10;
   int process_rank, process_num;
   // Initialize the MPI environment
   //&argc, &argv
   MPI_Init(NULL, NULL);
   // Find out rank, size
   
   MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
   
   MPI_Comm_size(MPI_COMM_WORLD, &process_num);

   // We are assuming 2 processes for this task
   if (process_num != 2) {
      fprintf(stderr, "World size must be two for %s\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   int local_mem_sum, global_mem_sum;
   int n_size =4;
   int arr[] = {5, 3, 4, 7};

   for (int i = 0; i < n_size; i++) {
        local_mem_sum += arr[i];
    }

   //address of the send buffer
   int ping_pong_count = 0;
   int stride = 0;
   //receiver
   int partner_rank = (process_rank + 1) % 2;
   // Can we have this MPI Collective Operation???
   MPI_Barrier(MPI_COMM_WORLD);
   double start,time;
   while (ping_pong_count < PING_PONG_LIMIT) {
      if (process_rank == ping_pong_count % 2) {
         
         start = MPI_Wtime();
         // Increment the ping pong count before you send it
         ping_pong_count++; //+2

         MPI_Send(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
         time += MPI_Wtime() - start;
         printf("process %d sent and incremented ping_pong_count %d\n",
               process_rank, ping_pong_count);
         
      } else {
         if(stride + process_rank < PING_PONG_LIMIT){
            int recv_sum;
            start = MPI_Wtime();
            MPI_Recv(&recv_sum, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            //take the recv value and add it to the local mem
            local_mem_sum += recv_sum;
            time += MPI_Wtime() - start;
            printf("process %d received ping_pong_count %d\n",
                  process_rank, ping_pong_count);
         }
         stride *= 2; //since the thread block takes two input elements into shared mem
      }
   }
   
   //global
   if (process_rank % 2 == 0) {
      global_mem_sum = local_mem_sum;
   }

   //

   //terminate MPI environment and releases all resources
   MPI_Finalize();
   printf("Time of sending the message back and forth is %d floating-point number of seconds\n", time / 10);
   printf("My rank is:%d\n",process_rank);
   printf("Array A sum of values is:%d\n",global_mem_sum);
   return 0;
}

#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include <chrono>


int main(int argc, char** argv){

    int cache_line_size;
    // in bytes
    if (argc > 1) cache_line_size = atoi(argv[1])/4;

    int num_threads;
    
    int size = 10e8;
    int i;
    int* a = (int*) malloc(size*sizeof(int));
    for (i = 0; i < size; i++) a[i] = 2;

    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    if (num_threads != 4) {
        printf("You are using %d threads. Please rerun with 4 threads!\n", num_threads);
        return 0;
    }

    int lres [4*cache_line_size];
    lres[0] = 0;
    lres[cache_line_size] = 0;
    lres[2*cache_line_size] = 0;
    lres[3*cache_line_size] = 0;

    auto start = std::chrono::steady_clock::now();

    int gres = 0;

    #pragma omp parallel private(i) shared(a,lres,size)
    {
        #pragma omp single 
        {
            #pragma omp task
            {
                for(i=0; i <size/4; i++)
                    lres[0] += a[i];
            }
            #pragma omp task
            {
                for(i=size/4; i <size/2; i++)
                    lres[cache_line_size] += a[i];
            }
            #pragma omp task
            {
                for(i=size/2; i <size/4*3; i++)
                    lres[2*cache_line_size] += a[i];
            }
            #pragma omp task
            {
                for(i=size/4*3; i <size; i++)
                    lres[3*cache_line_size] += a[i];
            }
            #pragma omp taskwait
            gres += lres[0];
            gres += lres[cache_line_size];
            gres += lres[2*cache_line_size];
            gres += lres[3*cache_line_size];
            printf("Parallel Solution = %d\n", gres);

        }
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time = std::chrono::duration<double>(0);
    time += end - start;
    printf("Time: %f\n", time.count());


    int seq_res = 0;
    start = std::chrono::steady_clock::now();

    for(i=0; i <size; i++)
        seq_res += a[i];

    printf("Sequential Solution = %d\n", seq_res);
    end = std::chrono::steady_clock::now();
    time = std::chrono::duration<double>(0);
    time += end - start;
    printf("Sequential Time: %f\n", time.count());

    return 0;
}
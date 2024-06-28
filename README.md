# group3_multicore_praktikal

## How to Run the code
- each folder has a different version of the code that shoud all be run similarly
- use `make` to compile the project source code
- to submit a job on superMUC use `sbatch job.scp`
- to check on current status of the job use `squeue -u <user_id>`
- to cancel a job use `scancel <job_id>`
- to check on output read `job.out`


##  Task 1-2
These two tasks share the same source code, the main difference is using hardware counters in the 2nd task with papi but that is mostly done in the command window.

## Task 3
This folder contains two version of source code, one version without changing the access pattern of the data in the relax jacobi algorithm, but focused on improving the speed of the residual calculations and data copying. in the Zick version it focues more on the access pattern but also with the fast swaping. These two versions both have advantages and disadvantages that are thoroughly explained in the report.

## Task 4 OpenMP Parallelization
This folder contains the source code used when parellelizing with OpenMP, it also covers the concept of the first touch. To run with different KMP_AFFINITY, one can add export commands in the job.scp.

## Task 5 MPI Parallelization
This folder contains two main folder, one using only MPI with blocking and non-blocking option, the other utilizes the idea of using a hyprid of MPI and OpenMP to get the most out of your hardware.
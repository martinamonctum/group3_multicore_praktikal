# Intel compiler
CC =  icc 
CFLAGS = -O3 -xhost -fno-alias
PAPI_INC = -I/dss/lrzsys/sys/spack/release/22.2.1/opt/skylake_avx512/papi/6.0.0.1-intel-pf2m6io/include
PAPI_LIB = -L/dss/lrzsys/sys/spack/release/22.2.1/opt/skylake_avx512/papi/6.0.0.1-intel-pf2m6io/lib
MPICC = mpicc

all: heat

heat : heat.o input.o misc.o timing.o relax_gauss.o relax_jacobi.o 
	$(CC) $(CFLAGS) ${PAPI_LIB} -o $@ $+ -lm -lpapi

%.o : %.c heat.h timing.h input.h
	$(CC) $(CFLAGS) $(PAPI_INC) -c -o $@ $< 

clean:
	rm -f *.o heat *~ *.ppm *.annot

remake : clean all

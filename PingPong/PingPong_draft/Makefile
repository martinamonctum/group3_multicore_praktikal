# Intel compiler
CC =  icc 
CFLAGS = -O0
MPICC = mpicc

all: main

main : main.o
	$(MPICC) $(CFLAGS) ${PAPI_LIB} -o $@ $+ -lm -lmpi

%.o : %.c 
	$(MPICC) $(CFLAGS) $(PAPI_INC) -c -o $@ $< 

clean:
	rm -f *.o *~ *.ppm *.annot *.optrpt

remake : clean all
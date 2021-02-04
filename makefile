CFLAGS = -g -Wall -Wstrict-prototypes
PROGS = make-matrix print-matrix mm-serial mm-parallel 
LDFLAGS = -lm
CC = gcc
MCC = mpicc

all: $(PROGS)

make-matrix: make-matrix.o matrix.o

print-matrix: print-matrix.o matrix.o

mm-serial: mm-serial.o matrix.o

mm-parallel: mm-parallel.c matrix.o graph_checkerboard_io.c
	$(MCC) $(LDFLAGS) -o $@ $+

clean: 
	rm -f $(PROGS) *.o core*
.POSIX:
.SUFFIX:
.SUFFIX: .c

CC = gcc
CFLAGS = -O2 -g
QCC = qcc
MPICC = mpicc
BASILISK = $(HOME)/basilisk/src
cylinder: _cylinder.c; $(CC) -o cylinder $(CFLAGS) -I$(BASILISK) _cylinder.c -lm
_cylinder.c: cylinder.c; $(QCC) -disable-dimensions $(QCCFLAGS) cylinder.c -source
deploy/cylinder.c: cylinder.c
	$(QCC) -disable-dimensions  -nolineno -source $(QCCFLAGS) cylinder.c && \
	mv _cylinder.c deploy/cylinder.c
deploy/cylinder_mpi.c: cylinder.c
	CC99=$(MPICC) $(QCC) -disable-dimensions  -nolineno -source $(QCCFLAGS) -D_MPI=1 cylinder.c && \
	mv _cylinder.c deploy/cylinder_mpi.c
dep: deploy/cylinder.c deploy/cylinder_mpi.c
clean:
	rm -f _cylinder.c cylinder

.POSIX:
.SUFFIX:
.SUFFIX: .c

CFLAGS = -O2 -g
QCC = qcc
MPICC = mpicc
BASILISK = $(HOME)/basilisk/src
cylinder: _cylinder.c; $(MPICC) -o cylinder $(CFLAGS) -I$(BASILISK) _cylinder.c -lm
_cylinder.c: cylinder.c; CC99=$(MPICC) $(QCC) -disable-dimensions -source $(QCCFLAGS) -D_MPI=1 cylinder.c
deploy/cylinder.c: cylinder.c
	CC99=$(MPICC) $(QCC) -disable-dimensions -source $(QCCFLAGS) -D_MPI=1 cylinder.c && \
	mv _cylinder.c deploy/cylinder_mpi.c
dep: deploy/cylinder.c deploy/cylinder_mpi.c
clean:
	rm -f _cylinder.c cylinder

_cylinder.c: output_htg.h  output_xdmf.h


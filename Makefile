.POSIX:
.SUFFIX:
.SUFFIX: .c

MPICC = mpicc
MPICCFLAGS = -O2 -g
QCC = qcc
all: _2 _3
_2.c: 2.c; $(QCC) $(QCCFLAGS) -D_MPI=1 2.c -source
_3.c: 3.c; $(QCC) $(QCCFLAGS) -D_MPI=1 3.c -source
.c:; $(MPICC) -o $@ $(MPICCFLAGS) $< -lm
clean:; rm -f _2 _2.c _3 _3.c

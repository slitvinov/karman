.POSIX:
.SUFFIX:
.SUFFIX: .c

MPICC = mpicc
MPICCFLAGS = -O2 -g
QCC = qcc
all: _main
_main.c: main.c; $(QCC) $(QCCFLAGS) -D_MPI=1 main.c -source
.c:; $(MPICC) -o $@ $(MPICCFLAGS) $< -lm
clean:; rm -f _main _main.c

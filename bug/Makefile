.POSIX:
.SUFFIXES:
.SUFFIXES: .c

CFLAGS = -O2 -g
QCC = qcc
MPICC = mpicc
BASILISK = $(HOME)/basilisk/src
main: _main.c; $(MPICC) -o main $(CFLAGS) -I$(BASILISK) _main.c -lm
_main.c: main.c; CC99=$(MPICC) $(QCC) -disable-dimensions -source $(QCCFLAGS) -D_MPI=1 main.c
clean:
	rm -f _main.c main

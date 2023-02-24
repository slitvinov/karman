.POSIX:
.SUFFIX:
.SUFFIX: .c

MPICC = mpicc -std=c99
QCC =  qcc
QCCFLAGS = -O2 -g
all: main
.c:; CC99='$(MPICC)' $(QCC) $(QCCFLAGS) -D_MPI=1 -o $@ $< -lm
clean:; rm -f main

.POSIX:
.SUFFIX:
.SUFFIX: .c

MPICC = mpicc.openmpi
QCC =  qcc
QCCFLAGS = -O2 -g
all: main
.c:; CC99=$(MPICC) $(QCC) $(QCCFLAGS) -D_MPI=1 -o $@ $< -lm

install:
	mkdir -p $(HOME)/.local/bin
	cp ppm2mp4 $(HOME)/.local/bin
clean:; rm -f main

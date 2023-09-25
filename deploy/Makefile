.POSIX:
.SUFFIXES:

CC = gcc
MPICC = mpicc
MPICCFLAGS = -O2 -g
CFLAGS = -O2
BASILISK = $(HOME)/basilisk
M = \
cylinder\
cylinder_mpi \

all: $M
cylinder_mpi: cylinder_mpi.c
	$(MPICC) $< -o $@ $(MPICCFLAGS) $(LDFLAGS) -I$(BASILISK) -lm
cylinder: cylinder.c
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS) -I$(BASILISK) -lm

clean:; rm -f -- $M
install: $B
	mkdir -p $(HOME)/.local/bin && \
	cp $B $(HOME)/.local/bin/


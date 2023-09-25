.POSIX:
.SUFFIXES:

CC = gcc
MPICC = mpicc
MPICCFLAGS = -O2 -g
CFLAGS = -O2
M = \
cylinder\
cylinder_mpi \

all: $M
cylinder_mpi: cylinder_mpi.c
	$(MPICC) $< -o $@ $(MPICCFLAGS) $(LDFLAGS) -lm
cylinder: cylinder.c
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS) -lm

clean:; rm -f -- $M
install: $B
	mkdir -p $(HOME)/.local/bin && \
	cp $B $(HOME)/.local/bin/


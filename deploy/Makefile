.POSIX:
.SUFFIXES:
.SUFFIXES: .c

CFLAGS = -g -O2
CFLAGS_OPENMP = -fopenmp
M = \
cylinder\

all: $M
.c:; $(CC) $(CFLAGS) $(CFLAGS_OPENMP) -o $@ $< $(LDFLAGS) -lm
clean:; rm -f -- $M

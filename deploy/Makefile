.POSIX:
.SUFFIXES:
.SUFFIXES: .c

CFLAGS = -O2
CFLAGS_OPENMP = -fopenmp
M = \
cylinder\

all: $M
.c:; $(CC) $< -o $@ $(CFLAGS) $(CFLAGS_OPENMP) $(LDFLAGS) -lm
clean:; rm -f -- $M

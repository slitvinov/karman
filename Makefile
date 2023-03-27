.POSIX:
.SUFFIX:
.SUFFIX: .c

BASILISK = $(HOME)/basilisk/src
CC = c99
MPICC = mpicc
MPICCFLAGS = -O2 -g
CFLAGS = -O2 -g
QCC = qcc
all: 2 3 sphere cylinder
_2.c: 2.c; $(QCC) $(QCCFLAGS) -D_MPI=1 2.c -source
_3.c: 3.c; $(QCC) $(QCCFLAGS) -D_MPI=1 3.c -source
_sphere.c: sphere.c; $(QCC) $(QCCFLAGS) -D_MPI=1 sphere.c -source
_cylinder.c: cylinder.c; $(QCC) $(QCCFLAGS) cylinder.c -source

2: _2.c; $(MPICC) -o $@ $(MPICCFLAGS) _2.c -lm
3: _3.c; $(MPICC) -o $@ $(MPICCFLAGS) _3.c $(BASILISK)/gl/libglutils.a $(BASILISK)/gl/libfb_osmesa.a  `pkg-config --libs osmesa glu` -lm
sphere: _sphere.c; $(MPICC) -o sphere $(MPICCFLAGS) _sphere.c $(BASILISK)/gl/libglutils.a $(BASILISK)/gl/libfb_osmesa.a  `pkg-config --libs osmesa glu` -lm
cylinder: _cylinder.c; $(CC) -o cylinder $(CFLAGS) _cylinder.c -lm

deploy/cylinder.c: cylinder.c
	$(QCC) -nolineno -source $(QCCFLAGS) cylinder.c && mv _cylinder.c deploy/cylinder.c

clean:; rm -f _2.c _3.c _cylinder.c 2 3 cylinder deploy/cylinder.c

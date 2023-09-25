.POSIX:
.SUFFIX:
.SUFFIX: .c

CC = c99
CFLAGS = -O2 -g
QCC = qcc
cylinder: _cylinder.c; $(CC) -o cylinder $(CFLAGS) _cylinder.c -lm
_cylinder.c: cylinder.c; $(QCC) $(QCCFLAGS) cylinder.c -source
deploy/cylinder.c: cylinder.c
	$(QCC) -nolineno -source $(QCCFLAGS) cylinder.c && \
	mv _cylinder.c deploy/cylinder.c

clean:
	rm -f _cylinder.c

.POSIX:
.SUFFIXES:
.SUFFIXES: .c
PREFIX = $(HOME)/.local
CFLAGS = -O3 -g
M = \
dump_2iso \
dump2xdmf \
dump_info \
dump_select \
stl2dump \
morton \

all: $M
.c:; $(CC) $< $(CFLAGS) -o $@
stl2dump: stl2dump.c
	$(CC) $< -fopenmp $(CFLAGS) -o $@ -lm
clean:; -rm $M

install: all
	mkdir -p -- '$(PREFIX)/bin' && cp -- sdf2lmp split '$(PREFIX)/bin'

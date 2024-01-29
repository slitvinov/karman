.POSIX:
.SUFFIXES:
.SUFFIXES: .c
PREFIX = $(HOME)/.local
CFLAGS = -O3 -g
M = \
dump_info \
dump_select \

all: $M
.c:; $(CC) $< $(CFLAGS) $(LDFLAGS) -o $@
clean:; -rm $M

install: all
	mkdir -p -- '$(PREFIX)/bin' && cp -- sdf2lmp split '$(PREFIX)/bin'

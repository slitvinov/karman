#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FREAD(ptr, size, nmemb)                                                \
  if (fread(ptr, size, nmemb, input_file) != (uint64_t)(nmemb)) {              \
    fprintf(stderr, "dump_info: error: fail to read from '%s'\n", input_path); \
    exit(1);                                                                   \
  }

#define FSEEK(offset)                                                          \
  if (fseek(input_file, offset, SEEK_CUR) != 0) {                              \
    fprintf(stderr, "dump_info: error: fseek failed in '%s'\n", input_path);   \
    exit(1);                                                                   \
  }

static int is_leaf(unsigned flags) {
  enum { leaf = 1 << 1 };
  return (leaf & leaf) != 0;
}

int main(int argc, char **argv) {
  int Verbose;
  long i, j;
  unsigned len;
  char *input_path, name[1024];
  FILE *input_file;
  struct coord {
    double x, y, z;
  };
  struct DumpHeader {
    double t;
    long len;
    int i, depth, npe, version;
    struct coord n;
  } header;

  Verbose = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr,
              "Usage: dump_info [-h] [-v] file.dump\n"
              "Options:\n"
              "  -h                          Print help message and exit\n"
              "  -v                          Verbose\n"
              "  file.dump                   basilisk dump\n");
      exit(1);
    case 'v':
      Verbose = 1;
      break;
    default:
      fprintf(stderr, "dump_info: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  if ((input_path = argv[0]) == NULL) {
    fprintf(stderr, "dump_info: error: file.dump xois not given\n");
    exit(1);
  }

  if ((input_file = fopen(input_path, "r")) == NULL) {
    fprintf(stderr, "dump_info: error: fail to open '%s'\n", input_path);
    exit(1);
  }

  FREAD(&header, sizeof header, 1);
  fprintf(stderr, "verbose: %d\n", header.version);
  fprintf(stderr, "t: %g\n", header.t);
  fprintf(stderr, "len: %ld\n", header.len);
  fprintf(stderr, "npe: %d\n", header.npe);
  fprintf(stderr, "depth: %d\n", header.depth);
  fprintf(stderr, "i: %d\n", header.i);
  fprintf(stderr, "n: [%g %g %g]\n", header.n.x, header.n.y, header.n.z);
  for (i = 0; i < header.len; i++) {
    FREAD(&len, sizeof len, 1);
    FREAD(&name, sizeof *name, len);
    name[len] = '\0';
    fprintf(stderr, "name[%d]: %s\n", len, name);
  }

  double o[4];
  unsigned flags, prev;
  double val, size;
  long offset, icell;
  int *pos;
  int level;

  FREAD(o, sizeof o, 1);
  fprintf(stderr, "origin: [%g %g %g]\n", o[0], o[1], o[2]);
  fprintf(stderr, "size: %g\n", o[3]);
  offset = header.len * sizeof(double) + sizeof(unsigned);
  if ((pos = calloc(1 + header.depth, sizeof *pos)) == NULL) {
    fprintf(stderr, "dump_info: error: calloc failed\n");
    exit(1);
  }

  level = 0;
  icell = 0;
  for (;;) {
    if (fread(&flags, sizeof flags, 1, input_file) != 1)
      break;
    for (i = 0; i < header.len; i++) {
      FREAD(&val, sizeof val, 1);
      if (i == 0)
        size = val;
    }
    if (icell > 0) {
      if (!is_leaf(flags) && is_leaf(prev))
        level--;
      else if (!is_leaf(flags) && !is_leaf(prev))
        level++;
    }
    fprintf(stderr, "tree size: %d %d %ld\n", level, is_leaf(flags),
            (long)size);
    prev = flags;
    icell++;
  }

  free(pos);
  if (fclose(input_file) != 0) {
    fprintf(stderr, "dump_info: error: fail to close '%s'\n", input_path);
    exit(1);
  }
}

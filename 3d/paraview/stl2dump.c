#include <assert.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Config {
  double R[3], L;
  int minlevel, maxlevel;
  char *stl_path, *dump_path;
};
struct Set {
  size_t M;
  int64_t *nodes;
};
struct coord {
  double x, y, z;
};
struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  struct coord n;
};
static long morton(uint32_t, uint32_t, uint32_t);
static int set_ini(size_t, void *, struct Set *);
static int set_add(struct Set *, int64_t);
static int set_has(struct Set *, int64_t);

enum { TABLE_DOUBLE, TABLE_INT, TABLE_PCHAR };
static struct {
  const char *name;
  int type;
  long offset;
} Table[] = {
    {"X0", TABLE_DOUBLE, offsetof(struct Config, R[0])},
    {"Y0", TABLE_DOUBLE, offsetof(struct Config, R[1])},
    {"Z0", TABLE_DOUBLE, offsetof(struct Config, R[2])},
    {"L", TABLE_DOUBLE, offsetof(struct Config, L)},
    {"minlevel", TABLE_INT, offsetof(struct Config, minlevel)},
    {"maxlevel", TABLE_INT, offsetof(struct Config, maxlevel)},
    {"stl_path", TABLE_PCHAR, offsetof(struct Config, stl_path)},
    {"dump_path", TABLE_PCHAR, offsetof(struct Config, dump_path)},
};
static char *fields[] = {"size", "cs"};

int main(int argc, char **argv) {
  char *end;
  FILE *stl_file, *dump_file;
  float *stl_ver;
  double lo[3], hi[3], r;
  int Verbose, d, j, level, status;
  struct Config config;
  int32_t stl_nt, i, inv_delta, x, y, z, u, v, w, ilo[3], ihi[3];
  uint64_t code, ncells;
  struct Set **set;
  void **work;
  size_t nbytes;
  struct DumpHeader header;

  Verbose = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr, "Usage: stl2dump [-h] [-v] X0 Y0 Z0 L minlevel maxlevel "
                      "file.stl basilisk.dump\n"
                      "Options:\n"
                      "  -h          print help message and exit\n"
                      "  -v          verbose\n");
      exit(1);
    case 'v':
      Verbose = 1;
      break;
    case '-':
      argv++;
      goto positional;
    default:
      fprintf(stderr, "stl2dump: error: unknown option '%s'\n", *argv);
      exit(1);
    }
positional:
  for (i = 0; i < sizeof(Table) / sizeof(*Table); i++) {
    if (*argv == NULL) {
      fprintf(stderr, "stl2dump: missing '%s' option\n", Table[i].name);
      exit(1);
    }
    switch (Table[i].type) {
    case TABLE_DOUBLE:
      *(double *)((void *)&config + Table[i].offset) = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "stl2dump: error: '%s' is not a double\n", *argv);
        exit(1);
      }
      break;
    case TABLE_INT:
      *(int *)((void *)&config + Table[i].offset) = strtol(*argv, &end, 10);
      if (*end != '\0') {
        fprintf(stderr, "stl2dump: error: '%s' is not an integer\n", *argv);
        exit(1);
      }
      break;
    case TABLE_PCHAR:
      *(char **)((void *)&config + Table[i].offset) = *argv;
      break;
    }
    argv++;
  }
  if ((stl_file = fopen(config.stl_path, "r")) == NULL) {
    fprintf(stderr, "stl2dump: error: fail to open '%s'\n", config.stl_path);
    exit(1);
  }
  if (fseek(stl_file, 80, SEEK_SET) != 0) {
    fprintf(stderr, "stl2dump: error: fail to read '%s'\n", config.stl_path);
    exit(1);
  }
  if (fread(&stl_nt, sizeof(stl_nt), 1, stl_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to read '%s'\n", config.stl_path);
    exit(1);
  }
  if (Verbose)
    fprintf(stderr, "stl2dump: stl_nt: %d\n", stl_nt);

  if ((stl_ver = malloc(9 * stl_nt * sizeof *stl_ver)) == NULL) {
    fprintf(stderr, "stl2dump: error: malloc failed\n");
    exit(1);
  }
  for (i = 0; i < stl_nt; i++) {
    fseek(stl_file, 3 * sizeof *stl_ver, SEEK_CUR);
    if (fread(&stl_ver[9 * i], sizeof *stl_ver, 9, stl_file) != 9) {
      fprintf(stderr, "stl2dump: error: fail to read '%s'\n", config.stl_path);
      exit(1);
    }
    fseek(stl_file, 2, SEEK_CUR);
  }
  if (fclose(stl_file) != 0) {
    fprintf(stderr, "stl2dump: error: fail to close '%s'\n", config.stl_path);
    exit(1);
  }

  set = malloc(config.maxlevel * sizeof *set);
  work = malloc(config.maxlevel * sizeof *work);
  nbytes = 1 << 28;
  for (i = 0; i < config.maxlevel; i++) {
    if ((work[i] = malloc(nbytes)) == NULL) {
      fprintf(stderr, "stl2dump: malloc failed\n");
      exit(1);
    }
    set[i] = malloc(sizeof *set[i]);
    set_ini(nbytes, work[i], set[i]);
  }

  ncells = 0;
  inv_delta = 1 << (config.maxlevel - 1);
  for (i = 0; i < stl_nt; i++) {
    fprintf(stderr, "i: %d\n", i);
    lo[0] = lo[1] = lo[2] = DBL_MAX;
    hi[0] = hi[1] = hi[2] = -DBL_MAX;
    for (j = 0; j < 3; j++) {
      for (d = 0; d < 3; d++) {
        r = stl_ver[9 * i + 3 * j + d];
        if (r < lo[d])
          lo[d] = r;
        if (r > hi[d])
          hi[d] = r;
      }
    }
    for (d = 0; d < 3; d++) {
      ilo[d] = (lo[d] - config.R[d]) / config.L * inv_delta;
      ihi[d] = (hi[d] - config.R[d]) / config.L * inv_delta + 2;
      if (ilo[d] < 0)
        ilo[d] = 0;
      if (ilo[d] > inv_delta)
        ilo[d] = inv_delta;
      if (ihi[d] < 0)
        ihi[d] = 0;
      if (ihi[d] > inv_delta)
        ihi[d] = inv_delta;
    }
    for (x = ilo[0]; x < ihi[0]; x++)
      for (y = ilo[1]; y < ihi[1]; y++)
        for (z = ilo[2]; z < ihi[2]; z++) {
          level = config.maxlevel - 1;
          u = x;
          v = y;
          w = z;
          for (;;) {
            code = morton(u, v, w);
            if (status = set_has(set[level], code)) {
              if (status == 2) {
                fprintf(stderr, "stl2dump: set problem: level: %d\n", level);
                exit(1);
              }
              break;
            }
            if (set_add(set[level], code) != 0) {
              fprintf(stderr, "stl2dump: set overflow: level: %d\n", level);
              exit(1);
            }
            // fprintf(stderr, "%ld\n", code);
            ncells++;
            if (level <= 0)
              break;
            u >>= 1;
            v >>= 1;
            w >>= 1;
            level--;
          }
        }
  }

  fprintf(stderr, "stl2dump: ncells: %ld\n", ncells);

  if ((dump_file = fopen(config.dump_path, "w")) == NULL) {
    fprintf(stderr, "stl2dump: error: fail to open '%s'\n", config.dump_path);
    exit(1);
  }

  header.t = 0;
  header.len = sizeof fields / sizeof *fields; /* */
  header.i = 0;
  header.depth = 7;
  header.npe = 1;
  header.version = 170901;
  header.n.x = 0;
  header.n.y = 0;
  header.n.z = 0;

  if (fwrite(&header, sizeof(header), 1, dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config.dump_path);
    exit(1);
  }
  unsigned len;
  for (i = 0; i < header.len; i++) {
    len = strlen(fields[i]);
    if (fwrite(&len, sizeof(len), 1, dump_file) != 1) {
      fprintf(stderr, "stl2dump: error: fail to write '%s'\n",
              config.dump_path);
      exit(1);
    }
    if (fwrite(fields[i], len, 1, dump_file) != 1) {
      fprintf(stderr, "stl2dump: error: fail to write '%s'\n",
              config.dump_path);
      exit(1);
    }
  }
  if (fwrite(config.R, sizeof(config.R), 1, dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config.dump_path);
    exit(1);
  }
  if (fwrite(&config.L, sizeof(config.L), 1, dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config.dump_path);
    exit(1);
  }
  //  traverse(0, 0, 0, l);

  if (fclose(dump_file) != 0) {
    fprintf(stderr, "stl2dump: error: fail to close '%s'\n", config.dump_path);
    exit(1);
  }
}

static long left(long x) {
  x = (x | x << 32) & 0x1f00000000ffffull;
  x = (x | x << 16) & 0x1f0000ff0000ffull;
  x = (x | x << 8) & 0x100f00f00f00f00full;
  x = (x | x << 4) & 0x10c30c30c30c30c3ull;
  x = (x | x << 2) & 0x1249249249249249ull;
  return x;
}

static long morton(uint32_t x, uint32_t y, uint32_t z) {
  return (left(z) << 2) | (left(y) << 1) | (left(x) << 0);
}

static int set_ini(size_t nbytes, void *memory, struct Set *set) {
  size_t i;
  set->M = nbytes / sizeof(int64_t);
  set->nodes = memory;
  for (i = 0; i < set->M; i++)
    set->nodes[i] = -1;
  return 0;
}
static int set_add(struct Set *set, int64_t key) {
  int x;
  size_t cnt;
  int64_t key0;
  assert(key >= 0);
  x = key % set->M;
  for (cnt = 0; cnt < set->M; cnt++) {
    key0 = set->nodes[x];
    if (key0 == -1) {
      set->nodes[x] = key;
      return 0;
    } else if (key0 == key) {
      return 0;
    }
    x = (x + 1) % set->M;
  }
  return 1;
}
static int set_has(struct Set *set, int64_t key) {
  int x;
  int64_t key0;
  size_t cnt;
  assert(key >= 0);
  x = key % set->M;
  for (cnt = 0; cnt < set->M; cnt++) {
    key0 = set->nodes[x];
    if (key0 == key) {
      return 1;
    } else if (key0 == -1) {
      return 0;
    }
    x = (x + 1) % set->M;
  }
  return 2;
}

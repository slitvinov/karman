#include <assert.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Set {
  size_t M;
  int64_t *nodes;
};
struct Config {
  double R[3], L;
  int minlevel, maxlevel;
  char *stl_path, *dump_path;
  FILE *dump_file;
  struct Set **set;
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
static uint64_t morton(uint64_t, uint64_t, uint64_t);
static int set_ini(size_t, void *, struct Set *);
static int set_add(struct Set *, int64_t);
static int set_has(struct Set *, int64_t);
static uint64_t traverse(uint64_t, uint64_t, uint64_t, int, struct Config *);

enum { TABLE_DOUBLE, TABLE_INT, TABLE_PCHAR };
static const struct {
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
static const int shift[][3] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
                               {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};
static char *fields[] = {"size", "cs", "level"};

int main(int argc, char **argv) {
  char *end;
  double lo[3], hi[3], r;
  FILE *stl_file;
  float *stl_ver;
  int64_t inv_delta;
  int32_t stl_nt, i, ilo[3], ihi[3];
  int Verbose, d, j, level;
  size_t nbytes;
  struct Config config;
  struct DumpHeader header;
  uint64_t code, ncells, u, v, w, x, y, z, size;
  unsigned len;
  void **work;

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

  config.set = malloc(config.maxlevel * sizeof *config.set);
  work = malloc(config.maxlevel * sizeof *work);
  nbytes = (1ul << 25) * sizeof *(*config.set)->nodes;
  for (i = 0; i < config.maxlevel; i++) {
    if ((work[i] = malloc(nbytes)) == NULL) {
      fprintf(stderr, "stl2dump: malloc failed\n");
      exit(1);
    }
    if ((config.set[i] = malloc(sizeof *config.set[i])) == NULL) {
      fprintf(stderr, "stl2dump: malloc failed\n");
      exit(1);
    }
    set_ini(nbytes, work[i], config.set[i]);
  }

  ncells = 0;
  inv_delta = 1ul << (config.maxlevel - 1);
  for (i = 0; i < stl_nt; i++) {
    if (i % 10 == 0)
      fprintf(stderr, "stl2dump: %d/%d\n", i, stl_nt);
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
            if (set_has(config.set[level], code))
              break;
            if (set_add(config.set[level], code) != 0) {
              fprintf(stderr, "stl2dump: set overflow: level: %d\n", level);
              exit(1);
            }
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
  if ((config.dump_file = fopen(config.dump_path, "w")) == NULL) {
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

  if (fwrite(&header, sizeof(header), 1, config.dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config.dump_path);
    exit(1);
  }
  for (i = 0; i < header.len; i++) {
    len = strlen(fields[i]);
    if (fwrite(&len, sizeof(len), 1, config.dump_file) != 1) {
      fprintf(stderr, "stl2dump: error: fail to write '%s'\n",
              config.dump_path);
      exit(1);
    }
    if (fwrite(fields[i], len, 1, config.dump_file) != 1) {
      fprintf(stderr, "stl2dump: error: fail to write '%s'\n",
              config.dump_path);
      exit(1);
    }
  }
  if (fwrite(config.R, sizeof(config.R), 1, config.dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config.dump_path);
    exit(1);
  }
  if (fwrite(&config.L, sizeof(config.L), 1, config.dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config.dump_path);
    exit(1);
  }
  size = traverse(0, 0, 0, 0, &config);
  fprintf(stderr, "stl2dump: size: %" PRIu64 "\n", size);

  for (i = 0; i < config.maxlevel; i++) {
    free(config.set[i]);
    free(work[i]);
  }
  free(work);
  free(config.set);
  free(stl_ver);
  if (fclose(config.dump_file) != 0) {
    fprintf(stderr, "stl2dump: error: fail to close '%s'\n", config.dump_path);
    exit(1);
  }
}

static uint64_t left(uint64_t x) {
  x = (x | x << 32) & 0x1f00000000ffffull;
  x = (x | x << 16) & 0x1f0000ff0000ffull;
  x = (x | x << 8) & 0x100f00f00f00f00full;
  x = (x | x << 4) & 0x10c30c30c30c30c3ull;
  x = (x | x << 2) & 0x1249249249249249ull;
  return x;
}

static uint64_t morton(uint64_t x, uint64_t y, uint64_t z) {
  return (left(z) << 2) | (left(y) << 1) | (left(x) << 0);
}

static int set_ini(size_t nbytes, void *memory, struct Set *set) {
  size_t i;
  set->M = nbytes / sizeof *set->nodes;
  set->nodes = memory;
  for (i = 0; i < set->M; i++)
    set->nodes[i] = -1;
  return 0;
}
static int set_add(struct Set *set, int64_t key) {
  int64_t key0;
  size_t cnt;
  uint64_t x;
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
    x = (x + 1 + cnt) % set->M;
  }
  fprintf(stderr, "stl2dump: error: set_add: over capacity\n");
  exit(1);
}
static int set_has(struct Set *set, int64_t key) {
  int64_t key0;
  size_t cnt;
  uint64_t x;
  if (key < 0) {
    fprintf(stderr, "stl2dump: error: key < 0\n");
    exit(1);
  }
  x = key % set->M;
  for (cnt = 0; cnt < set->M; cnt++) {
    key0 = set->nodes[x];
    if (key0 == key) {
      return 1;
    } else if (key0 == -1) {
      return 0;
    }
    x = (x + 1 + cnt) % set->M;
  }
  fprintf(stderr, "stl2dump: error: set_has failed\n");
  exit(1);
}

static uint64_t traverse(uint64_t x, uint64_t y, uint64_t z, int level,
                         struct Config *config) {
  double values[sizeof fields / sizeof *fields];
  int leaf, i;
  uint32_t leaf_code;
  uint64_t cell_size, u, v, w;
  long pos, curr, code;

  values[0] = 0;
  values[1] = 42;
  values[2] = level;
  code = morton(x, y, z);
  leaf = level >= config->minlevel &&
         (level == config->maxlevel || !set_has(config->set[level], code));
  leaf_code = leaf ? 2 : 0;
  if (fwrite(&leaf_code, sizeof(leaf_code), 1, config->dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config->dump_path);
    exit(1);
  }
  pos = ftell(config->dump_file);
  if (fwrite(values, sizeof(values), 1, config->dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config->dump_path);
    exit(1);
  }
  cell_size = 1;
  if (!leaf) {
    for (i = 0; i < sizeof shift / sizeof *shift; i++) {
      u = (x << 1) + shift[i][0];
      v = (y << 1) + shift[i][1];
      w = (z << 1) + shift[i][2];
      cell_size += traverse(u, v, w, level + 1, config);
    }
  }
  curr = ftell(config->dump_file);
  fseek(config->dump_file, pos, SEEK_SET);
  values[0] = cell_size;
  if (fwrite(&values[0], sizeof(values[0]), 1, config->dump_file) != 1) {
    fprintf(stderr, "stl2dump: error: fail to write '%s'\n", config->dump_path);
    exit(1);
  }
  fseek(config->dump_file, curr, SEEK_SET);
  return cell_size;
}

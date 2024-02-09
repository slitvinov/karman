#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FREAD(ptr, size, nmemb)                                                \
  if (fread(ptr, size, nmemb, context->input_file) != (uint64_t)(nmemb)) {     \
    fprintf(stderr, "dump_2iso: error: fail to read from '%s'\n",              \
            context->input_path);                                              \
    exit(1);                                                                   \
  }

#define FREAD0(ptr, size, nmemb)                                               \
  if (fread(ptr, size, nmemb, context.input_file) != (uint64_t)(nmemb)) {      \
    fprintf(stderr, "dump_2iso: error: fail to read from '%s'\n",              \
            context.input_path);                                               \
    exit(1);                                                                   \
  }

struct coord {
  double x, y, z;
};
struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  struct coord n;
};
static long traverse(int, void (*)(int, void *), void *);
static void process(int, void *);
static void counter(int, void *);
struct Context {
  struct DumpHeader header;
  int *index, malloc_level, m_level;
  double *values, X0, Y0, Z0, L0;
  FILE *input_file, *cells_file, *scalars_file, *field_file;
  char *input_path, *cells_path, *scalars_path, *field_path;
};
static const int shift[8][3] = {
    {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
    {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1},
};

int main(int argc, char **argv) {
  long i, input_offset;
  unsigned len;
  int Verbose;
  double o[4];
  char **names;
  struct Context context;
  Verbose = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(
          stderr,
          "Usage: dump_2iso [-h] [-v] file.dump in.cells in.scalar in.field\n"
          "Options:\n"
          "  -h                          Print help message and exit\n"
          "  -v                          Verbose\n"
          "  file.dump                   basilisk dump\n");
      exit(1);
    case 'v':
      Verbose = 1;
      break;
    default:
      fprintf(stderr, "dump_2iso: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  if ((context.input_path = argv[0]) == NULL) {
    fprintf(stderr, "dump_2iso: error: file.dump is not given\n");
    exit(1);
  }
  if ((context.cells_path = argv[1]) == NULL) {
    fprintf(stderr, "dump_2iso: error: in.cells is not given\n");
    exit(1);
  }
  if ((context.scalars_path = argv[2]) == NULL) {
    fprintf(stderr, "dump_2iso: error: in.scalars is not given\n");
    exit(1);
  }
  if ((context.field_path = argv[3]) == NULL) {
    fprintf(stderr, "dump_2iso: error: in.field is not given\n");
    exit(1);
  }
  if ((context.input_file = fopen(context.input_path, "r")) == NULL) {
    fprintf(stderr, "dump_2iso: error: fail to open '%s'\n",
            context.input_path);
    exit(1);
  }
  FREAD0(&context.header, sizeof context.header, 1);
  fprintf(stderr, "version: %d\n", context.header.version);
  fprintf(stderr, "t: %.16e\n", context.header.t);
  fprintf(stderr, "len: %ld\n", context.header.len);
  fprintf(stderr, "npe: %d\n", context.header.npe);
  fprintf(stderr, "depth: %d\n", context.header.depth);
  fprintf(stderr, "i: %d\n", context.header.i);
  fprintf(stderr, "n: [%g %g %g]\n", context.header.n.x, context.header.n.y,
          context.header.n.z);
  if ((names = malloc(context.header.len * sizeof *names)) == NULL) {
    fprintf(stderr, "dump_2iso: error: malloc failed\n");
    exit(1);
  }
  for (i = 0; i < context.header.len; i++) {
    FREAD0(&len, sizeof len, 1);
    names[i] = malloc((len + 1) * sizeof *names[i]);
    FREAD0(names[i], sizeof *names[i], len);
    names[i][len] = '\0';
    fprintf(stderr, "name: %s\n", names[i]);
  }
  FREAD0(o, sizeof o, 1);
  fprintf(stderr, "origin: [%g %g %g]\n", o[0], o[1], o[2]);
  context.X0 = o[0];
  context.Y0 = o[1];
  context.Z0 = o[2];
  context.L0 = o[3];
  fprintf(stderr, "size: %g\n", o[3]);
  context.malloc_level = 0;
  context.index = NULL;
  if ((context.values = malloc(context.header.len * sizeof *context.values)) ==
      NULL) {
    fprintf(stderr, "dump_2iso: error: malloc failed\n");
    exit(1);
  }
  context.m_level = 0;

  input_offset = ftell(context.input_file);
  traverse(0, counter, &context);
  fprintf(stderr, "m_level: %d\n", context.m_level);

  fseek(context.input_file, input_offset, SEEK_SET);
  if ((context.cells_file = fopen(context.cells_path, "w")) == NULL) {
    fprintf(stderr, "dump_2iso: error: fail to open '%s'\n",
            context.cells_path);
    exit(1);
  }
  if ((context.scalars_file = fopen(context.scalars_path, "w")) == NULL) {
    fprintf(stderr, "dump_2iso: error: fail to open '%s'\n",
            context.scalars_path);
    exit(1);
  }
  if ((context.field_file = fopen(context.field_path, "w")) == NULL) {
    fprintf(stderr, "dump_2iso: error: fail to open '%s'\n",
            context.field_path);
    exit(1);
  }
  traverse(0, process, &context);

  free(context.index);
  free(context.values);
  for (i = 0; i < context.header.len; i++)
    free(names[i]);
  free(names);
  if (fclose(context.input_file) != 0) {
    fprintf(stderr, "dump_2iso: error: fail to close '%s'\n",
            context.input_path);
    exit(1);
  }
  if (fclose(context.cells_file) != 0) {
    fprintf(stderr, "dump_2iso: error: fail to close '%s'\n",
            context.cells_path);
    exit(1);
  }
  if (fclose(context.scalars_file) != 0) {
    fprintf(stderr, "dump_2iso: error: fail to close '%s'\n",
            context.scalars_path);
    exit(1);
  }
  if (fclose(context.field_file) != 0) {
    fprintf(stderr, "dump_2iso: error: fail to close '%s'\n",
            context.field_path);
    exit(1);
  }
}

static void counter(int level, void *context_v) {
  struct Context *context;
  context = context_v;
  if (level > context->m_level)
    context->m_level = level;
}

static void process(int level, void *context_v) {
  int i;
  int Delta;
  int cell[4];
  float val;
  struct Context *context;
  context = context_v;
  cell[0] = 0;
  cell[1] = 0;
  cell[2] = 0;
  for (i = 1; i <= level; i++) {
    Delta = 1 << (context->m_level - i);
    cell[0] += Delta * shift[context->index[i]][0];
    cell[1] += Delta * shift[context->index[i]][1];
    cell[2] += Delta * shift[context->index[i]][2];
  }
  cell[3] = context->m_level - level;
  if (fwrite(cell, sizeof cell, 1, context->cells_file) != 1) {
    fprintf(stderr, "dump_2iso: error: fail to write '%s'\n",
            context->cells_path);
    exit(1);
  }

  val = context->values[8];
  if (fwrite(&val, sizeof(val), 1, context->scalars_file) != 1) {
    fprintf(stderr, "dump_2iso: error: fail to write '%s'\n",
            context->scalars_path);
    exit(1);
  }

  val = context->values[9];
  if (fwrite(&val, sizeof(val), 1, context->field_file) != 1) {
    fprintf(stderr, "dump_2iso: error: fail to write '%s'\n",
            context->field_path);
    exit(1);
  }
}

static long traverse(int level, void (*process)(int, void *), void *context_v) {
  enum { leaf = 2 };
  unsigned flags;
  long size, size0;
  struct Context *context;
  context = context_v;
  size0 = 1;
  FREAD(&flags, sizeof flags, 1);
  FREAD(context->values, sizeof *context->values, context->header.len);
  size = context->values[0];
  if (flags & leaf)
    process(level, context);
  if (flags & leaf) {
    /* */
  } else {
    while (level + 1 >= context->malloc_level) {
      context->malloc_level = 2 * context->malloc_level + 2;
      context->index = realloc(context->index,
                               context->malloc_level * sizeof *context->index);
    }
    for (context->index[level + 1] = 0; context->index[level + 1] < 8;
         context->index[level + 1]++)
      size0 += traverse(level + 1, process, context);
  }
  assert(size0 == size);
  return size;
}

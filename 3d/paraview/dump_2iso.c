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
  long nleaf;
  FILE *input_file;
  char *input_path;
};

static const double shift[8][3] = {
    {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0},
    {1, 0, 0}, {1, 0, 1}, {1, 1, 1}, {1, 1, 0},
};

int main(int argc, char **argv) {
  long i, input_offset;
  unsigned len;
  char *input_path;
  int Verbose;
  double o[4];
  char **names;
  struct Context context;
  Verbose = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr,
              "Usage: dump_2iso [-h] [-v] file.dump\n"
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
  if ((input_path = argv[0]) == NULL) {
    fprintf(stderr, "dump_2iso: error: file.dump xois not given\n");
    exit(1);
  }

  if ((context.input_file = fopen(input_path, "r")) == NULL) {
    fprintf(stderr, "dump_2iso: error: fail to open '%s'\n", input_path);
    exit(1);
  }
  FREAD0(&context.header, sizeof context.header, 1);
  fprintf(stderr, "verbose: %d\n", context.header.version);
  fprintf(stderr, "t: %g\n", context.header.t);
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

  context.nleaf = 0;
  fseek(context.input_file, input_offset, SEEK_SET);
  traverse(0, process, &context);

  fprintf(stderr, "nleaf: %ld\n", context.nleaf);
  free(context.index);
  free(context.values);
  for (i = 0; i < context.header.len; i++)
    free(names[i]);
  free(names);
  if (fclose(context.input_file) != 0) {
    fprintf(stderr, "dump_2iso: error: fail to close '%s'\n", input_path);
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
  int Delta, x, y, z;
  struct Context *context;
  context = context_v;

  x = 0;
  y = 0;
  z = 0;
  for (i = 1; i <= level; i++) {
    Delta = 1 << (context->m_level - i);
    x += Delta * shift[context->index[i]][0];
    y += Delta * shift[context->index[i]][1];
    z += Delta * shift[context->index[i]][2];
  }
  Delta = 1 << level;
  fprintf(stderr, "%d: %d: [%d %d %d]\n", level, Delta, x, y, z);
  context->nleaf++;
}

static long traverse(int level, void (*process)(int, void *), void *context_v) {
  enum { leaf = 2 };
  unsigned flags;
  long size, size0;
  struct Context *context;
  context = context_v;

  FREAD(&flags, sizeof flags, 1);
  FREAD(context->values, sizeof *context->values, context->header.len);
  size = context->values[0];
  size0 = 1;
  // if (flags & leaf)
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

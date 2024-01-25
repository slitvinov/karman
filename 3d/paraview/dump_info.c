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

static FILE *input_file;
static char *input_path;
struct coord {
  double x, y, z;
};
struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  struct coord n;
};
static struct DumpHeader header;
static int *index;
static int m_level;
static long traverse(int);
static double X0, Y0, Z0, L0;
static const double shift[8][3] = {
    {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0},
    {1, 0, 0}, {1, 0, 1}, {1, 1, 1}, {1, 1, 0},
};
static FILE *xyz_file;
static long ncell_total;

int main(int argc, char **argv) {
  FILE *file;
  long i, j;
  unsigned len;
  char *input_path, name[1024];
  int Verbose;
  double o[4];
  char xyz_path[FILENAME_MAX], attr_path[FILENAME_MAX], xdmf_path[FILENAME_MAX],
      *xyz_base, *attr_base;
  char *path = "a";
  snprintf(xyz_path, sizeof xyz_path, "%s.xyz.raw", path);
  snprintf(attr_path, sizeof attr_path, "%s.attr.raw", path);
  snprintf(xdmf_path, sizeof xdmf_path, "%s.xdmf2", path);
  xyz_base = xyz_path;
  attr_base = attr_path;
  for (j = 0; xyz_path[j] != '\0'; j++) {
    if (xyz_path[j] == '/' && xyz_path[j + 1] != '\0') {
      xyz_base = &xyz_path[j + 1];
      attr_base = &attr_path[j + 1];
    }
  }

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
  FREAD(o, sizeof o, 1);
  fprintf(stderr, "origin: [%g %g %g]\n", o[0], o[1], o[2]);
  X0 = o[0];
  Y0 = o[1];
  Z0 = o[2];
  L0 = o[3];
  fprintf(stderr, "size: %g\n", o[3]);
  m_level = 0;
  index = NULL;
  if ((xyz_file = fopen(xyz_path, "w")) == NULL) {
    fprintf(stderr, "%s:%d: fail to open '%s'\n", __FILE__, __LINE__, xyz_path);
    return 1;
  }
  ncell_total = 0;
  traverse(0);
  if (fclose(xyz_file) != 0) {
    fprintf(stderr, "dump_info: error: fail to close '%s'\n", xyz_path);
    exit(1);
  }
  free(index);
  if (fclose(input_file) != 0) {
    fprintf(stderr, "dump_info: error: fail to close '%s'\n", input_path);
    exit(1);
  }
  if ((file = fopen(xdmf_path, "w")) == NULL) {
    fprintf(stderr, "%s:%d: fail to open '%s'\n", __FILE__, __LINE__,
            xdmf_path);
    return 1;
  }
  fprintf(file,
          "<Xdmf\n"
          "    Version=\"2\">\n"
          "  <Domain>\n"
          "    <Grid>\n"
          "      <Topology\n"
          "          TopologyType=\"Hexahedron\"\n"
          "          Dimensions=\"%ld\"/>\n"
          "      <Geometry>\n"
          "        <DataItem\n"
          "            Dimensions=\"%ld 3\"\n"
          "            Format=\"Binary\">\n"
          "          %s\n"
          "        </DataItem>\n"
          "      </Geometry>\n",
          ncell_total, 8 * ncell_total, xyz_base);
  fprintf(file, "    </Grid>\n"
                "  </Domain>\n"
                "</Xdmf>\n");
  fclose(file);
}

static void process(int level) {
  int i, j;
  double Delta, x, y, z;
  float xyz[8 * 3];
  x = 0;
  y = 0;
  z = 0;
  assert(level < header.depth + 2);
  for (i = 1; i <= level; i++) {
    Delta = L0 * (1. / (1 << i));
    x += Delta * (shift[index[i]][0] - 0.5);
    y += Delta * (shift[index[i]][1] - 0.5);
    z += Delta * (shift[index[i]][2] - 0.5);
  }
  j = 0;
  for (i = 0; i < 8; i++) {
    xyz[j++] = x + Delta * (shift[i][0] - 0.5);
    xyz[j++] = y + Delta * (shift[i][1] - 0.5);
    xyz[j++] = z + Delta * (shift[i][2] - 0.5);
  }
  if (fwrite(xyz, sizeof xyz, 1, xyz_file) != 1) {
    fprintf(stderr, "dump_info: failed to write\n");
    exit(1);
  }
  // fprintf(stderr, ": %g %g %g: %d %d %d\n", x, y, z, index[1], index[2],
  // index[3]);
  fprintf(stdout, "%d %d %d\n", index[1], index[2], index[3]);
  ncell_total++;
}

static long traverse(int level) {
  enum {
    leaf = 1 << 1,
  };
  unsigned flags;
  double val;
  long size, size0;
  int i;
  FREAD(&flags, sizeof flags, 1);
  for (i = 0; i < header.len; i++) {
    FREAD(&val, sizeof val, 1);
    if (i == 0)
      size = val;
  }
  size0 = 1;
  if (level == 3)
    process(level);
  if (flags & leaf) {
    /* */
  } else {
    while (level + 1 >= m_level) {
      m_level = 2 * m_level + 1;
      index = realloc(index, m_level * sizeof *index);
    }
    for (index[level + 1] = 0; index[level + 1] < 8; index[level + 1]++)
      size0 += traverse(level + 1);
  }
  assert(size0 == size);
  return size;
}

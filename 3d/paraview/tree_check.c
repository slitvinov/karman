@include <stdint.h>
@include <stdlib.h>
@include <string.h>
#include "grid/octree.h"
#include "utils.h"
#include "output.h"
static scalar phi[];
int main(int argc, char **argv) {
  char *dump_path;
  FILE *dump_file;
  int Verbose;
  Verbose = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr, "Usage tree_check basilisk.dump\n");
      exit(1);
    case 'v':
      Verbose = 1;
      break;
    default:
      fprintf(stderr, "tree_check: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  if ((dump_path = *argv) == NULL) {
    fprintf(stderr, "tree_check: error: basilisk.dump is not given\n");
    exit(1);
  }
  if (Verbose)
    fprintf(stderr, "tree_check: starting on %d ranks\n", npe());

  if ((dump_file = fopen(dump_path, "r")) == NULL) {
    fprintf(stderr, "tree_check: error: failed to open '%s'\n", dump_path);
    exit(1);
  }
  init_grid(1 << 3);
  dump("dump");
  restore(fp = dump_file, list = {phi});
  if (fclose(dump_file) != 0) {
    fprintf(stderr, "tree_check: error: failed to close '%s'\n", dump_path);
    exit(1);
  }
  if (Verbose)
    fprintf(stderr, "tree_check: done\n");
}

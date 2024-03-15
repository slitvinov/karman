#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  char *dump_path;

  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr, "Usage tree_check basilisk.dump\n");
      exit(1);
    default:
      fprintf(stderr, "tree_check: error: unknown option '%s'\n", *argv);
      exit(1);
    }

  if ((dump_path = *argv) == NULL) {
    fprintf(stderr, "tree_check: error: basilisk.dump is not given\n");
    exit(1);
  }
  restore(dump_path);
}

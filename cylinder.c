@include <stdbool.h>
@include <stdint.h>
@include <string.h>
@include <stdint.h>
#include "embed.h"
#include "navier-stokes/centered.h"
#include "output_xdmf.h"
static const double diameter = 0.125;
static double reynolds, tend;
static int maxlevel, minlevel, period, Image, Surface, Verbose;
static const char *output_prefix;
u.n[left] = dirichlet(1);
p[left] = neumann(0);
pf[left] = neumann(0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);
u.n[embed] = fabs(y) > 0.45 ? neumann(0) : dirichlet(0);
u.t[embed] = fabs(y) > 0.45 ? neumann(0) : dirichlet(0);
face vector muv[];
int main(int argc, char **argv) {
  char *end;
  int MaxLevelFlag, MinLevelFlag, PeriodFlag, ReynoldsFlag, TendFlag;
  ReynoldsFlag = 0;
  MaxLevelFlag = 0;
  MinLevelFlag = 0;
  PeriodFlag = 0;
  Image = 0;
  Verbose = 0;
  TendFlag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(
          stderr,
          "Usage: cylinder [-h] [-i] [-v] -r <Reynolds "
          "number> -l <resolution level> -p <dump period> "
          "-e <end time>\n"
          "Options:\n"
          "  -h     Display this help message\n"
          "  -v     Verbose\n"
          "  -i     Enable PPM image dumping\n"
          "  -r <Reynolds number>     the Reynolds number (a decimal number)\n"
          "  -l <num>    Minimum resolution level (positive integer)\n"
          "  -m <num>    Maximum resolution level (positive integer)\n"
          "  -o <preifx>              a prefix for the output files\n"
          "  -p <dump period>         the dump period (positive integer)\n"
          "  -e <end time>            end time of the simulation (decimal "
          "number)\n"
          "\n"
          "Example usage:\n"
          "  ./cylinder -v -i -r 100 -l 7 -m 10 -p 100 -e 2 -o h\n");
      exit(1);
    case 'r':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error:  -r needs an argument\n");
        exit(1);
      }
      reynolds = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "cylinder: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      ReynoldsFlag = 1;
      break;
    case 'm':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -m needs an argument\n");
        exit(1);
      }
      maxlevel = strtol(*argv, &end, 10);
      if (*end != '\0' || maxlevel <= 0) {
        fprintf(stderr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      MaxLevelFlag = 1;
      break;
    case 'l':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -l needs an argument\n");
        exit(1);
      }
      minlevel = strtol(*argv, &end, 10);
      if (*end != '\0' || minlevel <= 0) {
        fprintf(stderr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      MinLevelFlag = 1;
      break;
    case 'p':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -p needs an argument\n");
        exit(1);
      }
      period = strtol(*argv, &end, 10);
      if (*end != '\0' || period <= 0) {
        fprintf(stderr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      PeriodFlag = 1;
      break;
    case 'i':
      Image = 1;
      break;
    case 'v':
      Verbose = 1;
      break;
    case 'e':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -e needs an argument\n");
        exit(1);
      }
      tend = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "cylinder: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      TendFlag = 1;
      break;
    case 'o':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -o needs an argument\n");
        exit(1);
      }
      output_prefix = *argv;
      break;
    default:
      fprintf(stderr, "cylinder: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  if (!ReynoldsFlag) {
    fprintf(stderr, "cylinder: error: -r is not set\n");
    exit(1);
  }
  if (!MinLevelFlag) {
    fprintf(stderr, "cylinder: error: -l is not set\n");
    exit(1);
  }
  if (!MaxLevelFlag) {
    fprintf(stderr, "cylinder: error: -m is not set\n");
    exit(1);
  }
  if (!PeriodFlag) {
    fprintf(stderr, "cylinder: error: -p is not set\n");
    exit(1);
  }
  if (!TendFlag) {
    fprintf(stderr, "cylinder: error: -e must be set\n");
    exit(1);
  }
  L0 = 4;
  origin(-0.5, -L0 / 2.);
  init_grid(1 << minlevel);
  mu = muv;
  run();
}
event properties(i++) { foreach_face() muv.x[] = fm.x[] * diameter / reynolds; }

event init(t = 0) {
  vertex scalar phi[];
  foreach_vertex() {
    double p0;
    p0 = 0.5 - y;
    p0 = min(p0, 0.5 + y);
    p0 = min(p0, sq(x) + sq(y) - sq(diameter / 2));
    phi[] = p0;
  }
  fractions(phi, cs, fs);
  foreach () {
    u.x[] = 0;
    u.y[] = 0;
  }
}

event dump(i++; t <= tend) {
  char path[FILENAME_MAX];
  scalar omega[], m[];
  FILE *fp;
  long nx, ny;
  coord n, b;

  if (i % period == 0) {
    if (Verbose) {
      fields_stats();
      if (pid() == 0)
        fprintf(stderr, "cylinder: %d: %09d %.16e\n", npe(), i, t);
    }
    if (output_prefix != NULL) {
      vorticity(u, omega);
      sprintf(path, "%s.%09d", output_prefix, i);
      output_xdmf(t, {p, omega}, {u}, path);
      if (Image) {
        foreach ()
          m[] = cs[] - 0.5;
        sprintf(path, "%s.%09d.png", output_prefix, i);
        output_ppm(omega, file = path, n = 2048,
                   box = {{-0.5, -0.5}, {L0 - 0.5, 0.5}}, min = -2 / diameter,
                   max = 2 / diameter, linear = false, mask = m);
      }
    }
  }
}
event adapt(i++) {
  adapt_wavelet({cs, u}, (double[]){1e-2, 3e-3, 3e-3}, maxlevel = maxlevel,
                minlevel = minlevel);
}

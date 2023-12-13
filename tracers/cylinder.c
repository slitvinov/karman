#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include "embed.h"
#include "navier-stokes/centered.h"
#include "output_xdmf.h"
#include "tracer.h"
static const double diameter = 0.5;
static const int minlevel = 7;
static double reynolds, tend;
static int maxlevel, period, Image, Surface, Verbose;
static const char *output_prefix;
scalar f[];
scalar * tracers = {f};
u.n[left] = dirichlet(1);
p[left] = neumann(0);
pf[left] = neumann(0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);
u.n[embed] = fabs(y) > 0.45 ? neumann(0) : dirichlet(0);
u.t[embed] = fabs(y) > 0.45 ? neumann(0) : dirichlet(0);
f[left]    = dirichlet(y < 0);
face vector muv[];
int main(int argc, char **argv) {
  char *end;
  int ReynoldsFlag;
  int LevelFlag;
  int PeriodFlag;
  int TendFlag;
  ReynoldsFlag = 0;
  LevelFlag = 0;
  PeriodFlag = 0;
  Image = 0;
  Verbose = 0;
  TendFlag = 0;
  output_prefix = NULL;
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
          "  -l <resolution level>    the resolution level (positive integer)\n"
          "  -p <dump period>         the dump period (positive integer)\n"
          "  -e <end time>            end time of the simulation (decimal "
          "  -o <preifx>              a prefix for the .raw output files\n"
          "number)\n"
          "\n"
          "Example usage:\n"
          "  ./cylinder -v -i -r 100 -l 10 -p 100 -e 2\n");
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
    case 'l':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -l needs an argument\n");
        exit(1);
      }
      maxlevel = strtol(*argv, &end, 10);
      if (*end != '\0' || maxlevel <= 0) {
        fprintf(stderr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      LevelFlag = 1;
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
  if (!LevelFlag) {
    fprintf(stderr, "cylinder: error: -l is not set\n");
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
    p0 = max(p0, sq(3*diameter / 8) - sq(x) - sq(y));
    double p1;
    p1 = (diameter / 16) - y;
    p1 = min(p1, (diameter / 16) + y);
    p0 = max(p0, p1);
    double p2;
    p2 = (diameter / 16) - x;
    p2 = min(p2, (diameter / 16) + x);
    p0 = max(p0, p2);
    phi[] = p0;
  }
  fractions(phi, cs, fs);
  foreach () {
    u.x[] = 0;
    u.y[] = 0;
  }
}

event dump(i++; t <= tend) {
  static long iframe = 0;
  char png[FILENAME_MAX], xdmf[FILENAME_MAX];
  scalar omega[], m[];
  FILE *fp;
  long nx, ny;
  coord n, b;

  if (iframe % period == 0) {
    if (Verbose) {
      fields_stats();
      if (pid() == 0)
        fprintf(stderr, "cylinder: %d: %09d %.16e\n", npe(), i, t);
    }

    vorticity(u, omega);
    if (output_prefix != NULL) {
      sprintf(xdmf, "%s.%09ld", output_prefix, iframe);
      vorticity(u, omega);
      output_xdmf({p, omega, f}, {u}, xdmf);
    }
    if (Image) {
      foreach ()
        m[] = cs[] - 0.5;
      sprintf(png, "%09ld.ppm", iframe);
      output_ppm(omega, file = png, n = 512,
                 box = {{-0.5, -0.5}, {L0 - 0.5, 0.5}}, min = -2 / diameter,
                 max = 2 / diameter, linear = false, mask = m);
      sprintf(png, "f.%09ld.ppm", iframe);
      output_ppm(f, file = png, n = 512,
                 box = {{-0.5, -0.5}, {L0 - 0.5, 0.5}}, min = 0,
                 max = 1,mask = m);
    }
  }
  iframe++;
}
event adapt(i++) {
  adapt_wavelet({cs, u, f}, (double[]){1e-2, 3e-3, 3e-3, 3e-2}, maxlevel = maxlevel,
                minlevel = minlevel);
}

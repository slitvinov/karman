#include <stdbool.h>
#include <stdint.h>
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "output_htg.h"
static const char *force_path;
static const double diameter = 0.2;
static const int minlevel = 7;
static double reynolds, tend;
static int maxlevel, period, Surface, Verbose;

u.n[left] = dirichlet(1);
p[left] = neumann(0);
pf[left] = neumann(0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

face vector muv[];
scalar vof[];
int main(int argc, char **argv) {
  char *end;
  int ReynoldsFlag, LevelFlag, PeriodFlag, TendFlag;
  ReynoldsFlag = 0;
  LevelFlag = 0;
  PeriodFlag = 0;
  Verbose = 0;
  TendFlag = 0;
  force_path = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(
          stderr,
          "Usage: cylinder [-h] [-i] [-v] [-f force file] -r <Reynolds "
          "number> -m <maximum resolution level> -p <dump period> "
          "-e <end time>\n"
          "Options:\n"
          "  -h     Display this help message\n"
          "  -v     Verbose\n"
          "  -r <Reynolds number>     the Reynolds number (a decimal number)\n"
          "  -m <resolution level>    the maximum resolution level (positive "
          "integer)\n"
          "  -p <dump period>         the dump period (positive integer)\n"
          "  -e <end time>            end time of the simulation (decimal "
          "number)\n"
          "  -f <force file>          force file\n"
          "\n"
          "Example usage:\n"
          "  ./cylinder -v -r 100 -m 10 -p 100 -e 2\n"
          "  ./cylinder -v -r 100 -m 10 -p 100 -e 2 -f force.dat\n");
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
    case 'f':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -f needs an argument\n");
        exit(1);
      }
      force_path = *argv;
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
  size(3.0);
  origin(-0.5, -L0 / 2.0, -L0 / 2.0);
  init_grid(1 << minlevel);
  mu = muv;
  run();
}
event properties(i++) { foreach_face() muv.x[] = fm.x[] * diameter / reynolds; }
event velocity(i++) {
  foreach ()
    foreach_dimension() u.x[] = vof[] * u.x[];
}
event init(t = 0) {
  vertex scalar phi[];
  foreach_vertex() {
    double p0;
    p0 = sq() + sq(y) - sq(diameter / 2);
    phi[] = p0;
  }
  refine(sq(x) + sq(y) < sq(1.05 * diameter / 2) &&
         sq(x) + sq(y) > sq(0.95 * diameter / 2) && level < maxlevel);
  fractions(phi, vof);
  foreach () {
    u.x[] = vof[];
    u.y[] = 0;
    u.z[] = 0;
  }
}

event dump(i++; t <= tend) {
  char png[FILENAME_MAX], htg[FILENAME_MAX];
  double fx, fy, fz;
  long nx, ny;
  scalar omega[], m[];
  static FILE *fp;
  static long iframe = 0;
  if (iframe % period == 0) {
    if (Verbose) {
      fields_stats();
      if (pid() == 0)
        fprintf(stderr, "cylinder: %d: %09d %.16e\n", npe(), i, t);
    }
    sprintf(htg, "h.%09ld.htg", iframe);
    vorticity(u, omega);
    output_htg({p, omega, vof}, {u}, htg);
    if (force_path) {
      fx = 0;
      fy = 0;
      fz = 0;
      foreach (reduction(+ : fx), reduction(+ : fy), reduction(+ : fz)) {
        double dv = (1 - vof[]) * dv();
        fx += u.x[] * dv;
        fy += u.y[] * dv;
        fz += u.z[] * dv;
      }
      fx /= dt;
      fy /= dt;
      fz /= dt;
      if (pid() == 0) {
        if (fp == NULL) {
          if ((fp = fopen(force_path, "w")) == NULL) {
            fprintf(stderr, "stl: error: fail to open '%s'\n", force_path);
            exit(1);
          }
        } else {
          if ((fp = fopen(force_path, "a")) == NULL) {
            fprintf(stderr, "stl: error: fail to open '%s'\n", force_path);
            exit(1);
          }
        }
        fprintf(fp, "%ld %.16e %.16e %.16e %.16e %.16e\n", iframe, t, fx,
                fy, fz, dt);
        fflush(fp);
      }
    }
  }
  iframe++;
}

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grid/octree.h"
#include "fractions.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "output_xdmf.h"
#include "predicate.h"
#include "predicate_c.h"
static const char *force_path, *output_prefix, *stl_path;
static const double diameter = 1;
static const int outlevel = 2;
static double reynolds, tend;
static int maxlevel, minlevel, period, Surface, Verbose;
static float *stl_ver;
static uint32_t stl_nt;
static int slice(double x, double y, double z, double Delta) {
  return z <= 0 && z + Delta >= 0;
}

static double vec_dot(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double edg_sq(const double a[3], const double b[3]) {
  double u[3];
  u[0] = a[0] - b[0];
  u[1] = a[1] - b[1];
  u[2] = a[2] - b[2];
  return vec_dot(u, u);
}

static double edg_point_distance2(const double a[3], const double b[3],
                                  const double p[3]) {
  enum { X, Y, Z };
  double t, s, x, y, z;

  s = edg_sq(a, b);
  if (s == 0)
    return edg_sq(p, a);
  t = ((b[X] - a[X]) * (p[X] - a[X]) + (b[Y] - a[Y]) * (p[Y] - a[Y]) +
       (b[Z] - a[Z]) * (p[Z] - a[Z])) /
      s;
  if (t > 1.0)
    return edg_sq(p, b);
  if (t < 0.0)
    return edg_sq(p, a);
  x = (1 - t) * a[X] + t * b[X] - p[X];
  y = (1 - t) * a[Y] + t * b[Y] - p[Y];
  z = (1 - t) * a[Z] + t * b[Z] - p[Z];
  return x * x + y * y + z * z;
}

static void vec_minus(const double a[3], const double b[3], /**/ double c[3]) {
  enum { X, Y, Z };
  c[X] = a[X] - b[X];
  c[Y] = a[Y] - b[Y];
  c[Z] = a[Z] - b[Z];
}

static double tri_point_distance2(const double a[3], const double b[3],
                                  const double c[3], const double p[3]) {
  enum { X, Y, Z };

  double u[3], v[3], q[3];
  double A, B, C, D, E, det;
  double t1, t2;
  double x, y, z;
  double d1, d2;

  vec_minus(b, a, u);
  vec_minus(c, a, v);
  B = vec_dot(v, u);
  E = vec_dot(u, u);
  C = vec_dot(v, v);
  det = B * B - E * C;
  if (det == 0) {
    d1 = edg_point_distance2(a, b, p);
    d2 = edg_point_distance2(b, c, p);
    if (d1 < d2)
      return d1;
    return d2;
  }
  vec_minus(a, p, q);
  A = vec_dot(v, q);
  D = vec_dot(u, q);
  t1 = (D * C - A * B) / det;
  t2 = (A * E - D * B) / det;
  if (t1 < 0)
    return edg_point_distance2(a, c, p);
  if (t2 < 0)
    return edg_point_distance2(a, b, p);
  if (t1 + t2 > 1)
    return edg_point_distance2(b, c, p);
  x = q[X] + t1 * u[X] + t2 * v[X];
  y = q[Y] + t1 * u[Y] + t2 * v[Y];
  z = q[Z] + t1 * u[Z] + t2 * v[Z];
  return x * x + y * y + z * z;
}

double embed_interpolate3(Point point, scalar s, coord p) {
  int i = sign(p.x), j = sign(p.y);
  int k = sign(p.z);
  x = fabs(p.x);
  y = fabs(p.y);
  z = fabs(p.z);
  if (cs[i] && cs[0, j] && cs[i, j] && cs[0, 0, k] && cs[i, 0, k] &&
      cs[0, j, k] && cs[i, j, k]) {
    return (((s[] * (1. - x) + s[i] * x) * (1. - y) +
             (s[0, j] * (1. - x) + s[i, j] * x) * y) *
                (1. - z) +
            ((s[0, 0, k] * (1. - x) + s[i, 0, k] * x) * (1. - y) +
             (s[0, j, k] * (1. - x) + s[i, j, k] * x) * y) *
                z);
  } else {
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i])
        val += fabs(p.x) * (s[i] - s[]);
      else if (cs[-i])
        val += fabs(p.x) * (s[] - s[-i]);
    }
    return val;
  }
}

trace void embed_force3(scalar p, vector u, face vector mu, coord *Fp,
                        coord *Fmu) {
  coord Fps = {0}, Fmus = {0};
  foreach (reduction(+ : Fps) reduction(+ : Fmus), nowarning)
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      double area = embed_geometry(point, &b, &n);
      area *= pow(Delta, dimension - 1);
      double Fn = area * embed_interpolate3(point, p, b);
      foreach_dimension() Fps.x += Fn * n.x;
    }
  *Fp = Fps;
  *Fmu = Fmus;
}

scalar f[];
face vector muv[];
scalar *tracers = {f};

u.n[left] = dirichlet(1);
p[left] = neumann(0);
pf[left] = neumann(0);
f[left] = dirichlet(y < 0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

int main(int argc, char **argv) {
  char *end;
  int ReynoldsFlag, MaxLevelFlag, MinLevelFlag, PeriodFlag, TendFlag;
  uint32_t stl_i;
  FILE *stl_file;
  ReynoldsFlag = 0;
  MaxLevelFlag = 0;
  MinLevelFlag = 0;
  PeriodFlag = 0;
  Verbose = 0;
  TendFlag = 0;
  output_prefix = NULL;
  force_path = NULL;
  stl_path = NULL;
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
          "  -l <resolution level>    the minimum resolution level (positive "
          "integer)\n"
          "  -m <resolution level>    the maximum resolution level (positive "
          "integer)\n"
          "  -o <preifx>              a prefix for the output files\n"
          "  -p <dump period>         the dump period (positive integer)\n"
          "  -e <end time>            end time of the simulation (decimal "
          "number)\n"
          "  -f <force file>          force file\n"
          "  -s <STL file>            geomtry file\n"
          "\n"
          "Example usage:\n"
          "  ./cylinder -v -r 100 -l 7 -m 10 -p 100 -e 2\n"
          "  ./cylinder -v -r 100 -l 7 -m 10 -p 100 -e 2 -f force.dat\n");
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
    case 's':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -s needs an argument\n");
        exit(1);
      }
      stl_path = *argv;
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
  if (!MaxLevelFlag) {
    fprintf(stderr, "cylinder: error: -m is not set\n");
    exit(1);
  }
  if (!MinLevelFlag) {
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
  if (stl_path) {
    if ((stl_file = fopen(stl_path, "r")) == NULL) {
      fprintf(stderr, "cylinder: error: fail to open '%s'\n", stl_path);
      exit(1);
    }
    fseek(stl_file, 80, SEEK_SET);
    if (fread(&stl_nt, sizeof(stl_nt), 1, stl_file) != 1) {
      fprintf(stderr, "cylinder: error: fail to doubled '%s'\n", stl_path);
      exit(1);
    }
    stl_ver = malloc(9 * stl_nt * sizeof *stl_ver);
    if (Verbose && pid() == 0)
      fprintf(stderr, "cylinder: triangles in STL file: %d\n", stl_nt);
    for (stl_i = 0; stl_i < stl_nt; stl_i++) {
      fseek(stl_file, 3 * sizeof *stl_ver, SEEK_CUR);
      if (fread(&stl_ver[9 * stl_i], sizeof *stl_ver, 9, stl_file) != 9) {
        fprintf(stderr, "cylinder: error: fail to read '%s'\n", stl_path);
        exit(1);
      }
      fseek(stl_file, 2, SEEK_CUR);
    }
    if (fclose(stl_file) != 0) {
      fprintf(stderr, "cylinder: error: fail to close '%s'\n", stl_path);
      exit(1);
    }
  }
  size(50);
  origin(-L0 / 2.5, -L0 / 2.0, -L0 / 2.0);
  init_grid(1 << outlevel);
  mu = muv;
  periodic(top);
  run();
}
event properties(i++) { foreach_face() muv.x[] = fm.x[] * diameter / reynolds; }
event init(t = 0) {
  refine(x < X0 + 0.8 * L0 && level < minlevel);
  if (stl_path) {
    vertex scalar phi[];
    phi.refine = phi.prolongation = fraction_refine;
    refine(sq(x) + sq(y) <= sq(diameter) && sq(x) + sq(y) >= sq(diameter / 2) &&
           level < maxlevel);
    predicate_ini();
    foreach_vertex() {
      if (sq(x) + sq(y) <= sq(1.25 * diameter / 2) &&
          sq(x) + sq(y) >= sq(0.75 * diameter / 2)) {
        double minimum;
        uint32_t intersect;
        uint32_t stl_i, j;
        double a[3], b[3], c[3], e[3], s[3], dist2;
        intersect = 0;
        minimum = DBL_MAX;
        for (stl_i = 0; stl_i < stl_nt; stl_i++) {
          j = 9 * stl_i;
          a[0] = stl_ver[j];
          a[1] = stl_ver[j + 1];
          a[2] = stl_ver[j + 2];

          b[0] = stl_ver[j + 3];
          b[1] = stl_ver[j + 4];
          b[2] = stl_ver[j + 5];

          c[0] = stl_ver[j + 6];
          c[1] = stl_ver[j + 7];
          c[2] = stl_ver[j + 8];

          s[0] = x;
          s[1] = y;
          s[2] = z;

          e[0] = s[0];
          e[1] = s[1];
          e[2] = s[2] + 2 * L0;

          dist2 = tri_point_distance2(a, b, c, s);
          if (dist2 < minimum)
            minimum = dist2;
          intersect += predicate_ray(s, e, a, b, c);
        }
        phi[] = intersect % 2 ? -dist2 : dist2;
      } else
        phi[] = 0.25 * diameter / 2;
    }
    if (Verbose && pid() == 0)
      fprintf(stderr, "cylinder: exported geomtry\n");
    free(stl_ver);
  } else {
    for (;;) {
      solid(cs, fs, sq(x) + sq(y) - sq(diameter / 2));
      astats s = adapt_wavelet({cs}, (double[]){3e-3}, maxlevel = maxlevel,
                               minlevel = minlevel);
      if (Verbose && pid() == 0)
        fprintf(stderr, "cylinder: refined %d cells\n", s.nf);
      if (s.nf == 0)
        break;
    }
    // fractions_cleanup(cs, fs);
  }
  foreach () {
    u.x[] = cs[];
    u.y[] = 0;
    u.z[] = 0;
  }
}
event velocity(i++; t <= tend) {
  char xdmf[FILENAME_MAX];
  coord Fp, Fmu;
  scalar omega[], m[];
  static FILE *fp;
  static long iframe = 0;
  if (iframe % period == 0) {
    if (Verbose) {
      fields_stats();
      if (pid() == 0)
        fprintf(stderr, "cylinder: %d: %09d %.16e %ld\n", npe(), i, t, grid->n);
    }
    if (output_prefix != NULL) {
      vorticity(u, omega);
      sprintf(xdmf, "%s.%09ld", output_prefix, iframe);
      output_xdmf({p, omega, f, cs}, {u}, NULL, xdmf);
      sprintf(xdmf, "%s.slice.%09ld", output_prefix, iframe);
      output_xdmf({p, omega, f, cs}, {u}, NULL, xdmf);
    }
    if (force_path) {
      embed_force3(p, u, mu, &Fp, &Fmu);
      if (pid() == 0) {
        if (fp == NULL) {
          if ((fp = fopen(force_path, "w")) == NULL) {
            fprintf(stderr, "cylinder: error: fail to open '%s'\n", force_path);
            exit(1);
          }
        } else {
          if ((fp = fopen(force_path, "a")) == NULL) {
            fprintf(stderr, "cylinder: error: fail to open '%s'\n", force_path);
            exit(1);
          }
        }
        fprintf(fp,
                "%ld %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e "
                "%.16e %.16e\n",
                iframe, t, Fp.x + Fmu.x, Fp.y + Fmu.y, Fp.z + Fmu.z, Fp.x, Fp.y,
                Fp.z, Fmu.x, Fmu.y, Fmu.z, dt);
        fflush(fp);
      }
    }
  }
  astats s = adapt_wavelet((scalar *){u}, (double[]){3e-2, 3e-2, 3e-2},
                           maxlevel = maxlevel, minlevel = minlevel);
  unrefine(!(x < X0 + 0.8 * L0));
  if (Verbose && iframe % period == 0 && pid() == 0)
    fprintf(stderr, "cylinder: refined %d cells, coarsened %d cells\n", s.nf,
            s.nc);
  iframe++;
}

#include "embed.h"
#include "navier-stokes/centered.h"
static const double d = 0.0625;
static int level;
static double reynolds;

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
u.n[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);
face vector muv[];

int main(int argc, char **argv) {
  char *end;
  int ReynoldsFlag;
  int LevelFlag;
  ReynoldsFlag = 0;
  LevelFlag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr, "usage: cylinder -r Reynolds\n");
      exit(1);
    case 'r':
      argv++;
      if (*argv == NULL) {
	fprintf(stderr, "cylinder: -r needs an argument\n");
	exit(1);
      }
      reynolds = strtod(*argv, &end);
      if (*end != '\0') {
	fprintf(stderr, "cylinder '%s' is not a number\n", *argv);
	exit(1);
      }
      ReynoldsFlag = 1;
      break;
    case 'l':
      argv++;
      if (*argv == NULL) {
	fprintf(stderr, "cylinder: -l needs an argument\n");
	exit(1);
      }
      level = strtol(*argv, &end, 10);
      if (*end != '\0' || level <= 0) {
	fprintf(stderr, "cylinder '%s' is not an integer\n", *argv);
	exit(1);
      }
      LevelFlag = 1;
      break;
    default:
      fprintf(stderr, "cylinder: unknown option '%s'\n", *argv);
      exit(1);
    }
  if (!ReynoldsFlag) {
    fprintf(stderr, "cylinder: -r is not set\n");
    exit(1);
  }
  if (!LevelFlag) {
    fprintf(stderr, "cylinder: -l is not set\n");
    exit(1);
  }
  L0 = 4;
  origin (-0.5, -L0/2.);
  N = 1024;
  mu = muv;
  run();
}
event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*d/reynolds;
}

event init (t = 0)
{
  vertex scalar phi[];
  foreach_vertex() {
    int i;
    double p0;
    p0 = 0.5 - y;
    p0 = min(p0, 0.5 + y);
    p0 = min(p0, sq(x) + sq(y) - sq(d/2));
    phi[] = p0;
  }
  fractions(phi, cs, fs);
}

event logfile (i += 10) {
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

event movies (i += 50; t <= 100)
{
  static long iframe = 0;
  char path[FILENAME_MAX];
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  sprintf(path, "%09ld.png", iframe);
  output_ppm (omega, file = path, box = {{-0.5, -0.45}, {3.5, 0.45}},
	      min = -200, max = 200, linear = false, mask = m);
  iframe++;
}
event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2, 3e-3, 3e-3}, level, 4);
}

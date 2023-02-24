#include "embed.h"
#include "navier-stokes/centered.h"
static const double d = 0.0625;
double Reynolds = 1100.;
int maxlevel = 14;
const double pos[][2] = {{0, 0}, {2.5, 5.0/6}, {5.0, 10.0/6}, {7.5, 15.0/6}};
int icase;

face vector muv[];
int main(int argc, char **argv) {
  L0 = 4;
  origin (-0.5, -L0/2.);
  N = 2048;
  mu = muv;
  if (!argv[1]) {
    fprintf(stderr, "karman1: needs an argument\n");
    exit(1);
  }
  icase = atoi(argv[1]);
  run();
}
event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*d/Reynolds;
}

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
u.n[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{
  vertex scalar phi[];
  foreach_vertex() {
    int i;
    double p0;
    p0 = 0.5 - y;
    p0 = min(p0, 0.5 + y);
    for (i = 0; i < icase; i++) {
      p0 = min(p0, sq(x - pos[i][0] * d) + sq(y - pos[i][1] * d) - sq(d/2));
      p0 = min(p0, sq(x - pos[i][0] * d) + sq(y + pos[i][1] * d) - sq(d/2));
    }
    phi[] = p0;
  }
  fractions(phi, cs, fs);
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}

event logfile (i += 10)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event movies (i += 50; t <= 100)
{
  char path[FILENAME_MAX];
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  sprintf(path, "%d.0.mp4", icase);
  output_ppm (omega, file = path, box = {{-d, -d}, {d, d}},
	      min = -200, max = 200, linear = false, mask = m);
  sprintf(path, "%d.1.mp4", icase);
  output_ppm (omega, file = path, box = {{-d, -3.5 * d}, {10 * d, 3.5 * d}},
	      min = -200, max = 200, linear = false, mask = m);
  sprintf(path, "%d.2.mp4", icase);
  output_ppm (omega, file = path, box = {{-0.5, -0.45}, {3.5, 0.45}},
	      min = -200, max = 200, linear = false, mask = m);
}
event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2, 3e-3, 3e-3, 3e-3}, maxlevel, 4);
}

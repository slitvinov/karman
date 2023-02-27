#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
static const double d = 0.0625;
static const int code[] = {1, 3, 5, 7};
static const double pos[][2] = {{0, 0}, {2.5, 5.0/6}, {5.0, 10.0/6}, {7.5, 15.0/6}};
static double Reynolds = 10;
static int maxlevel = 6;
static int icase;

u.n[left]  = dirichlet(1);
p[left]    = neumann(0);
pf[left]   = neumann(0);
u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
u.n[embed] = fabs(y) > 0.45 ? neumann(0) : dirichlet(0);
u.t[embed] = fabs(y) > 0.45 ? neumann(0) : dirichlet(0);
u.r[embed] = fabs(y) > 0.45 ? neumann(0) : dirichlet(0);
face vector muv[];

int main(int argc, char **argv) {
  if (!argv[1]) {
    fprintf(stderr, "main: needs an case number\n");
    exit(1);
  }
  icase = atoi(argv[1]);
  L0 = 1;
  origin (-0.5, -L0/2., -0.5);
  N = 64;
  mu = muv;
  run();
}
event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*d/Reynolds;
}

event init (t = 0)
{
  vertex scalar phi[];
  foreach_vertex() {
    int i;
    double p0;
    p0 = 0.5 - y;
    p0 = min(p0, 0.5 + y);
    p0 = min(p0, 0.125 - z);
    p0 = min(p0, 0.125 + z);
    for (i = 0; i < icase; i++) {
      p0 = min(p0, sq(x - pos[i][0] * d) + sq(y - pos[i][1] * d) - sq(d/2));
      p0 = min(p0, sq(x - pos[i][0] * d) + sq(y + pos[i][1] * d) - sq(d/2));
    }
    phi[] = p0;
  }
  fractions(phi, cs, fs);
}

event logfile (i += 10) {
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

event movies (i += 1; t <= 100)
{
  static long iframe = 0;
  char path[FILENAME_MAX];
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    omega[] = ((u.y[1] - u.y[-1]) - (u.x[0,1] - u.x[0,-1]))/(4.*Delta);
  view (fov = 15.65, quat = {-0.52,0.31,0.38,-0.7},
	tx = -0.045, ty = 0.015, width = 640, height = 480, bg = {0,0,0});  
  isosurface ("u.x", 0.5, color = "level");
  sprintf(path, "%d.2.%09ld.png", code[icase - 1], iframe);
  save(path);
  iframe++;
}

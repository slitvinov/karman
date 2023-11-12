#include <stdbool.h>
#include <stdint.h>
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "output_htg.h"
bool emerged = true;
scalar csm1[];
static double embed_interpolate3 (Point point, scalar s, coord b)
{
  int i = sign(b.x), j = sign(b.y);
  int k = sign(b.z);
  if (cs[i,0,0] && cs[0,j,0] && cs[i,j,0] &&
      cs[0,0,k] && cs[i,0,k] && cs[0,j,k] && cs[i,j,k] &&
      (emerged || (csm1[i,0,0] && csm1[0,j,0] && csm1[i,j,0] &&
		   csm1[0,0,k] && csm1[i,0,k] && csm1[0,j,k] && csm1[i,j,k]))) {
    double val_0, val_k;
    // bilinear interpolation in x-y-planes when all neighbors are defined
    val_0 = (s[0,0,0]*(1. - fabs(b.x)) + s[i,0,0]*fabs(b.x))*(1. - fabs(b.y)) +
      (s[0,j,0]*(1. - fabs(b.x)) + s[i,j,0]*fabs(b.x))*fabs(b.y);
    val_k = (s[0,0,k]*(1. - fabs(b.x)) + s[i,0,k]*fabs(b.x))*(1. - fabs(b.y)) +
      (s[0,j,k]*(1. - fabs(b.x)) + s[i,j,k]*fabs(b.x))*fabs(b.y);
    // trilinear interpolation when all neighbors are defined
    return (val_0*(1. - fabs(b.z)) + val_k*fabs(b.z));
  }
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(b.x);
      if (cs[i] &&
	  (emerged || (csm1[] && csm1[i])))
	val += fabs(b.x)*(s[i] - s[]);
      else if (cs[-i] &&
	       (emerged || (csm1[] && csm1[-i])))
	val += fabs(b.x)*(s[] - s[-i]);
    }
    return val;
  }
}

static
double embed_geometry3 (Point point, coord * b, coord * n)
{
  *n = facet_normal (point, cs, fs);
  double alpha = plane_alpha (cs[], *n);
  double area = plane_area_center (*n, alpha, b);
  normalize (n);
  return area;
}

static
coord embed_gradient3 (Point point, vector u, coord b, coord n)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet;
    double vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs, n, b, vb, &val);
      dudn.x += u.x[]*val; // For pathological situations
    }
    else // Neumann
      dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}

trace
void embed_force3 (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
  coord Fps = {0}, Fmus = {0};
  foreach (reduction(+:Fps) reduction(+:Fmus), nowarning)
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      double area = embed_geometry3 (point, &b, &n);
      area *= pow (Delta, dimension - 1);
      double Fn = area * embed_interpolate3 (point, p, b);
      foreach_dimension()
	Fps.x += Fn*n.x;
      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fs.x[] + fs.x[1];
	}
	mua /= (fa + SEPS);
	coord dudn = embed_gradient3 (point, u, b, n);
	foreach_dimension()
	  Fmus.x -= area*mua*(dudn.x*(sq (n.x) + 1.) +
			      dudn.y*n.x*n.y +
			      dudn.z*n.x*n.z);
      }
    }
  *Fp = Fps; *Fmu = Fmus;
}


static const char *force_path, *output_prefix;
static const double diameter = 2;
static double reynolds, tend;
static int maxlevel, minlevel, period, Surface, Verbose;

u.n[left] = dirichlet(1);
p[left] = neumann(0);
pf[left] = neumann(0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

face vector muv[];
int main(int argc, char **argv) {
  char *end;
  int ReynoldsFlag, MaxLevelFlag, MinLevelFlag, PeriodFlag, TendFlag;
  ReynoldsFlag = 0;
  MaxLevelFlag = 0;
  MinLevelFlag = 0;
  PeriodFlag = 0;
  Verbose = 0;
  TendFlag = 0;
  output_prefix = NULL;
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
	  "  -l <resolution level>    the minimum resolution level (positive "
	  "integer)\n"
	  "  -m <resolution level>    the maximum resolution level (positive "
	  "integer)\n"
	  "  -o <preifx>              a prefix for the output files\n"
	  "  -p <dump period>         the dump period (positive integer)\n"
	  "  -e <end time>            end time of the simulation (decimal "
	  "number)\n"
	  "  -f <force file>          force file\n"
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
  size(5);
  origin(-2, -L0 / 2.0, -L0 / 2.0);
  init_grid(1 << minlevel);
  mu = muv;
  run();
}
event properties(i++) { foreach_face() muv.x[] = fm.x[] * diameter / reynolds; }
event init(t = 0) {
  int l;
  double eps;
  vertex scalar phi[];
  for (l = minlevel + 1; l <= maxlevel; l++)
    refine(sq(x) + sq(y) <= sq(1.25 * diameter / 2) &&
	   sq(x) + sq(y) >= sq(0.95 * diameter / 2) &&
	   level < l);
  foreach_vertex() phi[] = sq(x) + sq(y) - sq(diameter / 2);
  fractions(phi, cs, fs);
  foreach () {
    u.x[] = cs[];
    u.y[] = 0;
    u.z[] = 0;
  }
}
event velocity(i++; t <= tend) {
  char htg[FILENAME_MAX];
  coord Fp, Fmu;
  long nx, ny;
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
      sprintf(htg, "%s.%09ld.htg", output_prefix, iframe);
      vorticity(u, omega);
      output_htg({p, omega, cs}, {u}, htg);
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
  astats s = adapt_wavelet((scalar*){u}, (double[]){3e-3, 3e-3, 3e-3},
			   maxlevel = maxlevel, minlevel = minlevel);
  if (Verbose && iframe % period == 0 && pid() == 0)
    fprintf(stderr, "cylinder: refined %d cells, coarsened %d cells\n", s.nf,
	    s.nc);
  iframe++;
}

event metric (i = 0)
{
  foreach() {
    cs[] = 1.;
    csm1[] = 1.;
  }
  foreach_face()
    fs.x[] = 1.;
#if TREE
  cs.restriction   = restriction_average;
  csm1.restriction = restriction_average;
  cs.refine        = embed_fraction_refine;
  cs.prolongation = fraction_refine;
  csm1.refine = csm1.prolongation = fraction_refine;
  foreach_dimension()
    fs.x.prolongation = embed_face_fraction_refine_x;
#endif
  restriction ({cs, csm1, fs});
  assert (is_constant (cm) || cm.i == cs.i);
  cm = cs;
  fm = fs;
  csm1.nodump = true;
}

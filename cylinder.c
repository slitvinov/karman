#include "embed.h"
#include "navier-stokes/centered.h"
static const double diameter = 0.125;
static double reynolds;
static int level, period, Image, Surface;
u.n[left] = dirichlet(1.);
p[left] = neumann(0.);
pf[left] = neumann(0.);
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
u.n[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);
face vector muv[];

int main(int argc, char **argv) {
  char *end;
  int ReynoldsFlag;
  int LevelFlag;
  int PeriodFlag;
  ReynoldsFlag = 0;
  LevelFlag = 0;
  PeriodFlag = 0;
  Image = 0;
  Surface = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr,
	      "usage: cylinder [-i] [-s] -r <Reynolds number> -l <resolution "
	      "level> -p <dump period>\n");
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
	fprintf(stderr, "cylinder: '%s' is not a positive integer\n", *argv);
	exit(1);
      }
      LevelFlag = 1;
      break;
    case 'p':
      argv++;
      if (*argv == NULL) {
	fprintf(stderr, "cylinder: -p needs an argument\n");
	exit(1);
      }
      period = strtol(*argv, &end, 10);
      if (*end != '\0' || period <= 0) {
	fprintf(stderr, "cylinder: '%s' is not a positive integer\n", *argv);
	exit(1);
      }
      PeriodFlag = 1;
      break;
    case 'i':
      argv++;
      Image = 1;
      break;
    case 's':
      argv++;
      Surface = 1;
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
  if (!PeriodFlag) {
    fprintf(stderr, "cylinder: -p is not set\n");
    exit(1);
  }
  L0 = 4;
  origin(-0.5, -L0 / 2.);
  N = 1024;
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

event dump(i++; t <= 100) {
  static long iframe = 0;
  char png[FILENAME_MAX], raw[FILENAME_MAX], xdmf[FILENAME_MAX], surface[FILENAME_MAX];
  scalar omega[], m[];
  FILE *fp;
  double sx, sy, xp, yp, ox, oy;
  float v;
  long k, j, nx, ny;
  char *names[] = {"ux", "uy", "p"};
  coord n, b;
  double omega_surface, theta;

  if (iframe % period == 0) {
    fprintf(stderr, "cylinder: %09d %.3e\n", i, t);
    vorticity(u, omega);
    sprintf(raw, "%09ld.raw", iframe);
    if ((fp = fopen(raw, "w")) == NULL) {
      fprintf(stderr, "cylinder: fail to write to '%s'\n", raw);
      exit(1);
    }
    ox = X0;
    oy = -0.5;
    nx = N;
    sx = L0 / nx;
    ny = 1.0 / sx;
    sy = 1.0 / ny;
    for (k = 0; k < ny; k++) {
      yp = oy + sy * k + sy / 2.;
      for (j = 0; j < nx; j++) {
	xp = ox + sx * j + sx / 2;
	v = interpolate(u.x, xp, yp);
	fwrite(&v, sizeof v, 1, fp);
	v = interpolate(u.y, xp, yp);
	fwrite(&v, sizeof v, 1, fp);
	v = interpolate(p, xp, yp);
	fwrite(&v, sizeof v, 1, fp);
      }
    }
    if (fclose(fp) != 0) {
      fprintf(stderr, "cylinder: fail to close '%s'\n", raw);
      exit(1);
    }

    sprintf(xdmf, "a.%09ld.xdmf2", iframe);
    if ((fp = fopen(xdmf, "w")) == NULL) {
      fprintf(stderr, "cylinder: fail to write to '%s'\n", xdmf);
      exit(1);
    }
    fprintf(fp, "\
<?xml version=\"1.0\" ?>\n\
<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n\
<Xdmf Version=\"2.0\">\n\
 <Domain>\n\
   <Grid Name=\"Grid\">\n\
     <Topology TopologyType=\"2DCORECTMesh\" Dimensions=\"%ld %ld\"/>\n\
     <Geometry GeometryType=\"ORIGIN_DXDY\">\n\
       <DataItem Name=\"Origin\" Dimensions=\"2\">%.16e %.16e</DataItem>\n\
       <DataItem Name=\"Spacing\" Dimensions=\"2\">%.16e %.16e</DataItem>\n\
     </Geometry>\n\
",
	    ny + 1, nx + 1, oy, ox, sy, sx);
    for (j = 0; j < sizeof names / sizeof *names; j++)
      fprintf(fp, "\
     <Attribute Name=\"%s\" Center=\"Cell\">\n\
	<DataItem ItemType=\"HyperSlab\" Dimensions=\"%ld %ld\">\n\
	  <DataItem Dimensions=\"3 2\">0 %ld 1 %ld %ld %ld</DataItem>\n\
	  <DataItem Dimensions=\"%ld %ld\" Format=\"Binary\">%s</DataItem>\n\
	</DataItem>\n\
     </Attribute>\n\
",
	      names[j], ny, nx, j, sizeof names / sizeof *names, ny, nx, ny, 3 * nx, raw);
    fprintf(fp,"\
   </Grid>\n\
 </Domain>\n\
</Xdmf>\n\
");
    if (fclose(fp) != 0) {
      fprintf(stderr, "cylinder: fail to close '%s'\n", xdmf);
      exit(1);
    }

    if (Surface) {
      sprintf(surface, "surface.%09ld.raw", iframe);
      if ((fp = fopen(surface, "w")) == NULL) {
	fprintf(stderr, "cylinder: fail to write to '%s'\n", surface);
	exit(1);
      }
      foreach(serial)
	if (cs[] > 0. && cs[] < 1.) {
	  embed_geometry(point, &b, &n);
	  omega_surface = embed_vorticity(point, u, b, n);
	  x += b.x*Delta, y += b.y*Delta;
	  theta = atan2(y, x);
	  fwrite(&theta, sizeof theta, 1, fp);
	  fwrite(&omega_surface, sizeof omega_surface, 1, fp);
	}
      if (fclose(fp) != 0) {
	fprintf(stderr, "cylinder: fail to close '%s'\n", surface);
	exit(1);
      }
    }
    if (Image) {
      foreach ()
	m[] = cs[] - 0.5;
      sprintf(png, "%09ld.ppm", iframe);
      output_ppm(omega, file = png, box = {{-0.5, -0.5}, {L0 - 0.5, 0.5}},
		 min = -5 / diameter, max = 5 / diameter, linear = false,
		 mask = m);
    }
  }
  iframe++;
}
event adapt(i++) {
  adapt_wavelet({cs, u}, (double[]){1e-2, 3e-3, 3e-3}, level, 4);
}

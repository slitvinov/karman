#include "embed.h"
#include "navier-stokes/centered.h"
static const double diameter = 0.125;
static double reynolds;
static int level, period, Image, Surface, Zoom, nzoom;
u.n[left] = dirichlet(1.);
p[left] = neumann(0.);
pf[left] = neumann(0.);
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
u.n[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);
face vector muv[];

int dump_fields(const char *raw, const char *xdmf, double t, double ox,
                double oy, double ex, double ey, long nx) {
  long k, j, ny;
  double sx, sy, xp, yp;
  float v;
  FILE *fp;
  char *names[] = {"ux", "uy", "p"};
  sx = ex / nx;
  ny = ey / sx;
  sy = ey / ny;
  if ((fp = fopen(raw, "w")) == NULL) {
    fprintf(stderr, "cylinder: fail to write to '%s'\n", raw);
    exit(1);
  }
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
    return 1;
  }
  if ((fp = fopen(xdmf, "w")) == NULL) {
    fprintf(stderr, "cylinder: fail to write to '%s'\n", xdmf);
    return 1;
  }
  fprintf(fp, "\
<?xml version=\"1.0\" ?>\n\
<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n\
<Xdmf Version=\"2.0\">\n\
 <Domain>\n\
   <Grid>\n\
     <Time Value=\"%.16e\"/>\n\
     <Topology TopologyType=\"2DCORECTMesh\" Dimensions=\"%ld %ld\"/>\n\
     <Geometry GeometryType=\"ORIGIN_DXDY\">\n\
       <DataItem Name=\"Origin\" Dimensions=\"2\">%.16e %.16e</DataItem>\n\
       <DataItem Name=\"Spacing\" Dimensions=\"2\">%.16e %.16e</DataItem>\n\
     </Geometry>\n\
",
          t, ny + 1, nx + 1, oy, ox, sy, sx);
  for (j = 0; j < sizeof names / sizeof *names; j++)
    fprintf(fp, "\
     <Attribute Name=\"%s\" Center=\"Cell\">\n\
	<DataItem ItemType=\"HyperSlab\" Dimensions=\"%ld %ld\">\n\
	  <DataItem Dimensions=\"3 2\">0 %ld 1 %ld %ld %ld</DataItem>\n\
	  <DataItem Dimensions=\"%ld %ld\" Format=\"Binary\">%s</DataItem>\n\
	</DataItem>\n\
     </Attribute>\n\
",
            names[j], ny, nx, j, sizeof names / sizeof *names, ny, nx, ny,
            3 * nx, raw);
  fprintf(fp, "\
   </Grid>\n\
 </Domain>\n\
</Xdmf>\n\
");
  if (fclose(fp) != 0) {
    fprintf(stderr, "cylinder: fail to close '%s'\n", xdmf);
    return 1;
  }
  return 0;
}

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
  Zoom = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr, "Usage: cylinder [-h] [-i] [-s] [-z <number of cells>] -r <Reynolds number> -l <resolution level> -p <dump period>\n"
	      "Options:\n"
	      "  -h     Display this help message\n"
	      "  -i     Enable PPM image dumping\n"
	      "  -s     Enable surface file dumping\n"
	      "  -z <number of cells>     Set the zoom level (positive integer)\n"
	      "  -r <Reynolds number>     Set the Reynolds number (a decimal number)\n"
	      "  -l <resolution level>    Set the resolution level (positive integer)\n"
	      "  -p <dump period>         Set the dump period (positive integer)\n"
	      "\n"
	      "Example usage:\n"
	      "  ./cylinder -i -s -r 100 -l 10 -p 100\n"
	      "  ./cylinder -z 1024 -r 100 -l 10 -p 100\n");
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
      level = strtol(*argv, &end, 10);
      if (*end != '\0' || level <= 0) {
        fprintf(stderr, "cylinder: error: '%s' is not a positive integer\n", *argv);
        exit(1);
      }
      LevelFlag = 1;
      break;
    case 'z':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -z needs an argument\n");
        exit(1);
      }
      nzoom = strtol(*argv, &end, 10);
      if (*end != '\0' || nzoom <= 0) {
        fprintf(stderr, "cylinder: error: '%s' is not a positive integer\n", *argv);
        exit(1);
      }
      Zoom = 1;
      break;
    case 'p':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "cylinder: error: -p needs an argument\n");
        exit(1);
      }
      period = strtol(*argv, &end, 10);
      if (*end != '\0' || period <= 0) {
        fprintf(stderr, "cylinder: error: '%s' is not a positive integer\n", *argv);
        exit(1);
      }
      PeriodFlag = 1;
      break;
    case 'i':
      Image = 1;
      break;
    case 's':
      Surface = 1;
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
  char png[FILENAME_MAX], raw[FILENAME_MAX], xdmf[FILENAME_MAX],
      surface[FILENAME_MAX];
  scalar omega[], m[];
  FILE *fp;
  long nx, ny;
  coord n, b;
  double omega_surface, theta;

  if (iframe % period == 0) {
    fprintf(stderr, "cylinder: %09d %.16e\n", i, t);

    sprintf(xdmf, "a.%09ld.xdmf2", iframe);
    sprintf(raw, "%09ld.raw", iframe);
    if (dump_fields(raw, xdmf, t, X0, -0.5, L0, 1.0, N) != 0)
      exit(1);

    if (Zoom) {
      sprintf(xdmf, "z.%09ld.xdmf2", iframe);
      sprintf(raw, "z.%09ld.raw", iframe);
      if (dump_fields(raw, xdmf, t, -1.25 * diameter, -1.25 * diameter,
                      2.5 * diameter, 2.5 * diameter, nzoom) != 0)
        exit(1);
    }

    if (Surface) {
      sprintf(surface, "surface.%09ld.raw", iframe);
      if ((fp = fopen(surface, "w")) == NULL) {
        fprintf(stderr, "cylinder: error: fail to write to '%s'\n", surface);
        exit(1);
      }
      foreach (serial)
        if (cs[] > 0. && cs[] < 1.) {
          embed_geometry(point, &b, &n);
          omega_surface = embed_vorticity(point, u, b, n);
          x += b.x * Delta, y += b.y * Delta;
          theta = atan2(y, x);
          fwrite(&theta, sizeof theta, 1, fp);
          fwrite(&omega_surface, sizeof omega_surface, 1, fp);
        }
      if (fclose(fp) != 0) {
        fprintf(stderr, "cylinder: error: fail to close '%s'\n", surface);
        exit(1);
      }
    }
    if (Image) {
      foreach ()
        m[] = cs[] - 0.5;
      sprintf(png, "%09ld.ppm", iframe);
      vorticity(u, omega);
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

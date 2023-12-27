#include "navier-stokes/centered.h"

int main()
{
  DT = 1;
  run();
}

event init (i = 0) {
  restore("restore");
}

event dumping (t = 2)
  dump();

event ending (t = 6);

event logfile (i++)
  fprintf (stderr, "%d %g %g\n", i, t, dt);

Compile

```
$ make
c99 -g -O2 -fopenmp -o cylinder cylinder.c  -lm
```

Run

```
$ ./cylinder -i -s -r 100 -l 10 -p 100
cylinder: 000000000 0.000e+00
cylinder: 000000100 1.721e-01
cylinder: 000000200 3.875e-01
cylinder: 000000300 6.143e-01
...
```

Plot velocity and pressure fields

```
$ python plot.py a.000000300.xdmf2
```

open `ux.png`, `uy.png`, or `p.png`.

Plot vorticity on the surface of the cylinder

```
$ python surface.py surface.000000400.raw
```

open `omega.png`

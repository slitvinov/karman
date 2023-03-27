Compile

```
$ c99 cylinder.c -O2 -o cylinder -lm
```

Run

```
$ ./cylinder -r 100 -l 10 -p 100 -i
cylinder: 000000000 0.000e+00
cylinder: 000000100 1.721e-01
cylinder: 000000200 3.875e-01
cylinder: 000000300 6.143e-01
...
```

Plot

```
$ python plot.py
```

open `ux.png` and `uy.png`.

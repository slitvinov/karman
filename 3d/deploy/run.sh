#!/bin/sh

# ../stl/cylinder.py
../stl/center.py $HOME/cfd.stl
 ../paraview/stl2dump -v -- -4.8 -6 -6 12    7 9   64 center.stl basilisk.dump
exec mpiexec ./cylinder -v -r 220 -l 7 -m 9 -p 10 -e 3000 -d basilisk.dump -o h

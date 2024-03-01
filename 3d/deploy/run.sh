#!/bin/sh

# ../stl/cylinder.py
(cd ../paraview && make stl2dump)
../stl/center.py $HOME/cfd.stl
../paraview/stl2dump -v -- -2.5 -3.125 -3.124 6.25  5 9  64 center.stl basilisk.dump
# exec mpiexec ./cylinder -v -r 220 -l 8 -m 11 -p 10 -e 3000 -d basilisk.dump -o h

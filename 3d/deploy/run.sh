#!/bin/sh

# ../stl/cylinder.py
(cd ../paraview && make stl2dump)
../stl/center.py $HOME/cfd.stl
../paraview/stl2dump -v -- -5 -6.25 -6.25 12.5  7 10  64 center.stl basilisk.dump
# exec mpiexec ./cylinder -v -r 220 -l 6 -m 10 -p 100 -e 3000 -d basilisk.dump -o h


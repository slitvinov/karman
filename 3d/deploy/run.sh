#!/bin/sh

# ../stl/cylinder.py
(cd ../paraview && make stl2dump)
../stl/center.py $HOME/cfd.stl
../paraview/stl2dump -v -- -20 -25 -25 50   8 11  64 center.stl basilisk.dump
exec mpiexec ./cylinder -v -r 220 -l 8 -m 11 -p 10 -e 3000 -d basilisk.dump -o h

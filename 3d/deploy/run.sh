#!/bin/sh

../stl/cylinder.py
../paraview/stl2dump -v -- -4.8 -6 -6 12    6 8   64 center.stl basilisk.dump
exec mpiexec ./cylinder -v -r 220 -l 6 -m 8 -p 10 -e 3000 -d basilisk.dump -o h

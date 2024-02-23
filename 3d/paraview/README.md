```
make
../stl/cylinder.py
./stl2dump -v -- -1 -1 -1 2   3 6      2 center.stl basilisk.dump
./dump2xdmf basilisk.dump basilisk
mpiexec.mpich -n 2 ../cylinder -v -r 100 -l 3 -m 6 -p 1 -e 3000 -f force.dat -z 2 -S cylinder -o h -d basilisk.dump
```

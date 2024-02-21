```
make
../stl/cylinder.py
./stl.py -v -- -1 -1 -1 2   3 6 center.stl basilisk.dump
./dump2xdmf basilisk.dump basilisk
../cylinder -v -r 100 -l 6 -m 8 -p 10 -e 3000 -f force.dat -z 10 -S cylinder -o h -d basilisk.dump
```

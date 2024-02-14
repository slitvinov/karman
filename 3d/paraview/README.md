```
make
../stl/cylinder.py
./stl2dump -v -- -3 -3 -3 6  4 7 center.stl  q.dump
./dump2xdmf -v q.dump q
paraview q.xdmf2 center.stl
```

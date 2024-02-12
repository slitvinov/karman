```
make
../stl/cylinder.py
./stl2dump -v -- -2 -2 -2 4  4 7 center.stl  q.dump
./dump2xdmf -v q.dump q
```

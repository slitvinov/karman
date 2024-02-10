import struct
import sys
import mmap
import itertools
import numpy as np


def read_stl(path):
    with open(path, "rb+") as file:
        mm = mmap.mmap(file.fileno(), 0)
        nt, = struct.unpack('<i', mm[80:80 + 4])
        return np.ndarray((nt, 3, 3), np.dtype("<f4"), mm, 80 + 4 + 12,
                          (36 + 12 + 2, 12, 4))


Verbose = 0
argv = iter(sys.argv)
while True:
    sys.argv.pop(0)
    if len(sys.argv) and len(sys.argv[0]) > 1 and sys.argv[0][0] == '-':
        if sys.argv[0][1] == 'h':
            sys.stderr.write(
                """usage: stl.py [-h] [-v] X0 Y0 Z0 L minlevel maxlevel file.stl basilisk.dump
Options:
  -h     display this help message
  -v     verbose
""")
            sys.exit(2)
        elif sys.argv[0][1] == 'v':
            Verbose = 1
        elif sys.argv[0][1] == '-':
            sys.argv.pop(0)
            break
        else:
            sys.stderr.write("stl.py: error: unknown option '%s'\n" %
                             sys.argv[0])
            sys.exit(2)
    else:
        break
try:
    X0 = float(sys.argv.pop(0))
    Y0 = float(sys.argv.pop(0))
    Z0 = float(sys.argv.pop(0))
    L = float(sys.argv.pop(0))
except IndexError:
    sys.stderr.write("stl.py: error: not enough arguments\n")
except ValueError as e:
    sys.stderr.write("stl.py: error: %s\n" % e)
    sys.exit(1)
try:
    minlevel = int(sys.argv.pop(0))
    maxlevel = int(sys.argv.pop(0))
    stl_path = sys.argv.pop(0)
    dump_path = sys.argv.pop(0)
except IndexError:
    sys.stderr.write("stl.py: error: not enough arguments\n")
    sys.exit(1)
except ValueError as e:
    sys.stderr.write("stl.py: error: %s\n" % e)
    sys.exit(1)

try:
    stl = read_stl(stl_path)
except FileNotFoundError as e:
    sys.stderr.write("stl.py: error: fail to open %s\n" % stl_path)
    sys.exit(1)
except ValueError:
    sys.stderr.write("stl.py: error: fail to read STL file %s\n" % stl_path)
    sys.exit(1)

R0 = X0, Y0, Z0
cells = set()
for tri in stl:
    a, b, c = tri
    hi = np.max(tri, 0)
    lo = np.min(tri, 0)
    lo = (lo - R0) / L
    hi = (hi - R0) / L
    inv_delta = 1 << (maxlevel - 1)
    lo = [max(0, int(r)) for r in lo * inv_delta]
    hi = [min(int(r) + 1, inv_delta) for r in hi * inv_delta]
    for cell in itertools.product(*map(range, lo, hi)):
        cells.add(cell)
print(len(cells))

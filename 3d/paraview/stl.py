import itertools
import mmap
import numpy as np
import os
import struct
import sys

shift = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1),
         (1, 1, 0), (1, 1, 1))


def read_stl(path):
    with open(path, "rb+") as file:
        mm = mmap.mmap(file.fileno(), 0)
        nt, = struct.unpack('<i', mm[80:80 + 4])
        return np.ndarray((nt, 3, 3), np.dtype("<f4"), mm, 80 + 4 + 12,
                          (36 + 12 + 2, 12, 4))


def traverse(r, level):
    values = [0] * nfields
    values[-1] = level
    leaf = level >= minlevel and (level == maxlevel
                                  or not (*r, level) in cells)
    fmt = "%dd" % nfields
    dump.write(struct.pack("I", 2 if leaf else 0))
    pos = dump.tell()
    dump.write(struct.pack(fmt, *values))
    cell_size = 1
    if not leaf:
        for s in shift:
            r0 = [(r << 1) + s for r, s in zip(r, s)]
            cell_size += traverse(r0, level + 1)
    curr = dump.tell()
    dump.seek(pos, os.SEEK_SET)
    dump.write(struct.pack("d", cell_size))
    dump.seek(curr, os.SEEK_SET)
    return cell_size


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
sys.stderr.write("stl.py: stl_nt: %ld\n" % len(stl))

R0 = X0, Y0, Z0
cells = set()
inv_delta = 1 << (maxlevel - 1)
for tri in stl:
    a, b, c = tri
    hi = np.max(tri, 0)
    lo = np.min(tri, 0)
    lo = (lo - R0) / L
    hi = (hi - R0) / L

    lo = [max(0, int(r)) for r in lo * inv_delta]
    lo = [min(inv_delta, r) for r in lo]

    hi = [min(int(r) + 2, inv_delta) for r in hi * inv_delta]
    hi = [max(0, r) for r in hi]
    for cell in itertools.product(*map(range, lo, hi)):
        cells.add((*cell, maxlevel - 1))
        level = maxlevel - 1
        while True:
            level -= 1
            if level < 0:
                break
            cell = tuple(cell >> 1 for cell in cell)
            if (*cell, level) in cells:
                break
            else:
                cells.add((*cell, level))
sys.stderr.write("stl.py: cells: %ld\n" % len(cells))
fields = "size", "cs", "level"
nfields = len(fields)
t = 0
i = 0
depth = 7  # TODO:
npe = 1
version = 170901
n = 0, 0, 0

with open(dump_path, "wb") as dump:
    dump.write(struct.pack("dl4i3d", t, nfields, i, depth, npe, version, *n))
    for field in fields:
        fmt = "I%ds" % len(field)
        dump.write(struct.pack(fmt, len(field), str.encode(field)))
    dump.write(struct.pack("4d", X0, Y0, Z0, L))
    sys.stderr.write("stl.py: cells: %ld\n" % traverse((0, 0, 0), 0))

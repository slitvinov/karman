import collections
import struct
import os


def indicator(x, y, z):
    return (x - x0)**2 / rx**2 + (y - y0)**2 / ry**2 + (z - z0)**2 / rz**2 - 1


def traverse(l):
    values = [0] * nfields
    leaf = l == maxlevel
    fmt = "%dd" % nfields
    dump.write(struct.pack("I", 2 if leaf else 0))
    pos = dump.tell()
    dump.write(struct.pack(fmt, *values))
    cell_size = 1
    if leaf:
        pass
    else:
        for index[l + 1] in range(0, 8):
            cell_size += traverse(l + 1)
    curr = dump.tell()
    dump.seek(pos, os.SEEK_SET)
    dump.write(struct.pack("d", cell_size))
    dump.seek(curr, os.SEEK_SET)
    return cell_size


shift = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1),
         (1, 1, 0), (1, 1, 1))
minlevel = 1
maxlevel = 5
x0 = 1 / 2
y0 = 1 / 2
z0 = 1 / 2
rx = 1 / 4
ry = 1 / 5
rz = 1 / 6
cs = 0
fields = "size", "cs"
nfields = len(fields)
t = 0
i = 0
depth = 7  # TODO:
npe = 1
version = 170901
n = 0, 0, 0
X0, Y0, Z0 = 0, 0, 0
L0 = 1
index = collections.defaultdict(int)
with open("gen.dump", "wb") as dump:
    dump.write(struct.pack("dl4i3d", t, nfields, i, depth, npe, version, *n))
    for field in fields:
        fmt = "I%ds" % len(field)
        dump.write(struct.pack(fmt, len(field), str.encode(field)))
    dump.write(struct.pack("4d", X0, Y0, Z0, L0))
    traverse(0)

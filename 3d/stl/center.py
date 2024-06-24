#!/usr/bin/env python3
import numpy as np
import mmap
import sys
import struct
import os


def create(path, nt):
    with open(path, "wb") as file:
        file.write(struct.pack('<80s', b"Generated by center.py"))
        file.write(struct.pack('<i', nt))
        file.seek(nt * (4 * 12 + 2) - 1, os.SEEK_CUR)
        file.write(b'\0')
    with open(path, "rb+") as file:
        mm = mmap.mmap(file.fileno(), 0)
        sys.stderr.write("center.py: %s\n" % path)
        return np.ndarray((nt, 3, 3), np.dtype("<f4"), mm, 80 + 4 + 12,
                          (36 + 12 + 2, 12, 4))


def read(path):
    with open(path, "rb") as file:
        mm = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
        nt, = struct.unpack('<i', mm[80:80 + 4])
        return np.ndarray((nt, 3, 3), np.dtype("<f4"), mm, 80 + 4 + 12,
                          (36 + 12 + 2, 12, 4))


r = read(sys.argv[1])
nt = len(r)
x = np.ravel(r[:, :, 0])
y = np.ravel(r[:, :, 1])
z = np.ravel(r[:, :, 2])

xlo = min(x)
xhi = max(x)
ylo = min(y)
yhi = max(y)
zlo = min(z)
zhi = max(z)

xc = (xlo + xhi) / 2
yc = (ylo + yhi) / 2
zc = (zlo + zhi) / 2

L = [xhi - xlo, yhi - ylo, zhi - zlo]
ix, iy, iz = sorted([0, 1, 2], key=lambda i: L[i])
scale = 1 / L[ix]
R = [(x - xc) * scale, (y - yc) * scale, (z - zc) * scale]
x, y, z = R[ix], R[iy], R[iz]
sys.stderr.write("center.py: %g %g\ncenter.py: %g %g\ncenter.py: %g %g\n" %
                 (min(x), max(x), min(y), max(y), min(z), max(z)))

r = create("center.stl", nt)
np.copyto(r[:, :, 0], np.reshape(x, (nt, -1)))
np.copyto(r[:, :, 1], np.reshape(y, (nt, -1)))
np.copyto(r[:, :, 2], np.reshape(z, (nt, -1)))

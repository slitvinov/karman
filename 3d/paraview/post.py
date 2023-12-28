import sys
import re
import numpy as np
import matplotlib.patches
import matplotlib.pyplot as plt

path = sys.argv[1]
dtype = np.dtype("float32")

path = re.sub("\.xdmf2$", "", path)
path = re.sub("\.attr\.raw$", "", path)
path = re.sub("\.xyz\.raw$", "", path)

xyz_path = path + ".xyz.raw"
attr_path = path + ".attr.raw"
xyz = np.memmap(xyz_path, dtype)
ncell = xyz.size // (3 * 8)
assert ncell * 3 * 8 == xyz.size
patches = []

attr = np.memmap(attr_path, dtype)
attr = attr.reshape((ncell, -1))
u = attr[:, 2:2 + 3]

for i in range(ncell):
    j = 8 * 3 * i
    x = xyz[j]
    y = xyz[j + 1]

    k = j + 3 * 7

    u = xyz[k] - x
    v = xyz[k + 1] - y
    patches.append(matplotlib.patches.Rectangle((x, y), u, v, fill=None))

L = 2.5
plt.axis((-L / 2 + L / 10, L / 2 + L / 10, -L / 2, L / 2))
plt.axis('scaled')
plt.gca().add_collection(
    matplotlib.collections.PatchCollection(patches, match_original=True))
plt.show()

import matplotlib.pylab as plt
import numpy as np
import os
import re
import sys


def read(path):
    with open(path) as file:
        for line in file:
            if re.match('^	  <DataItem Dimensions=.*Format="Binary"', line):
                raw = re.sub('^[ \t]*<[^<]*>', '', line)
                raw = re.sub('<.*\n$', '', raw)
                shape = re.sub('^	  <DataItem Dimensions="', '', line)
                shape = re.sub('".*\n$', '', shape)
                ny, nx = [int(shape) for shape in shape.split()]
                break
    dirname = os.path.dirname(path)
    with open(os.path.join(dirname, raw), "rb") as file:
        return np.ndarray((3, nx // 3, ny), 'float32', file.read(), order='F')


ux, uy, p = read(sys.argv[1])
for name, field in zip(("ux.png", "uy.png", "p.png"), (ux, uy, p)):
    plt.imshow(field.T)
    plt.savefig(name, bbox_inches='tight')

import matplotlib.pylab as plt
import numpy as np
import os
import re


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
		ux, uy = np.ndarray((2, nx // 2, ny),
		                    'float32',
		                    file.read(),
		                    order='F')
		return ux, uy


ux, uy = read('a.000000300.xdmf2')
plt.imshow(ux.T)
plt.savefig('ux.png', bbox_inches='tight')
plt.imshow(uy.T)
plt.savefig('uy.png', bbox_inches='tight')

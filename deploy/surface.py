import matplotlib.pylab as plt
import numpy as np
import sys


def read(path):
    nitems = 2
    dtype = np.dtype('float64')
    with open(path, 'rb') as file:
        buffer = file.read()
        n = len(buffer) // (dtype.itemsize * nitems)
        return np.ndarray((nitems, n), dtype, buffer, order='F')


theta, omega = read(sys.argv[1])
plt.plot(theta, omega, 'o')
plt.xlabel("theta, radians")
plt.ylabel("vorticity, 1 / second")
plt.savefig("omega.png")

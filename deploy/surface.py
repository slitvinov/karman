import matplotlib.pylab as plt
import numpy as np
import os
import re
import sys


def read(path):
    with open(path) as file:
        return np.ndarray((4, -1), 'float64', file.read(), order='F')


read(sys.argv[0])

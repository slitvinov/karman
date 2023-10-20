import numpy as np
import math

xp, yp = np.loadtxt("tritton.txt").T
xp = np.log(xp)
yp = np.log(yp)

x = np.linspace(np.min(xp), np.max(xp), 100)
y = np.interp(x, xp, yp)

for u, v in zip(x, y):
    u = math.exp(u)
    v = math.exp(v)
    print("%.16e %.16e" % (u, v))

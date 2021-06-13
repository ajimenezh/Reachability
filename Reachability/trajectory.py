import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

c1 = 160
c2 = -1.6
c3 = 53
c4 = -3.5
c5 = 156
c6 = 78
s = 15

def f(x, u):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]

    
    y = np.zeros(5)
    y[0] = s*math.cos(x3)
    y[1] = s*math.sin(x3)
    y[2] = x4
    y[3] = -c1/s*x4 - c2*x5 + c3*u
    y[4] = (-1 - c4/s/s)*x4 - c5/s*x5 + c6/s*u

    return y

x = np.asarray([3, -1.8, 0.0, 0.0, 0.0])

t = 3.2
delta_t = 0.05
n = int(t / delta_t)
xx = [x[0]]
yy = [x[1]]

u_max = 0.038

for i in range (n+10):
    u = u_max
    if i >= n/4 and i < 3*n/4:
        u = -u
    if i >= n:
        u = 0.0

    
    y = f(x, u)
    x = x + y*delta_t

    print (i, i*delta_t, x)

    xx.append(x[0])
    yy.append(x[1])

patches = []
rect = mpatches.Rectangle([27.5, -4], 5, 2, ec="none")
patches.append(rect)
#label(grid[1], "Rectangle")

fig, ax = plt.subplots()

plt.plot(xx, yy)

plt.ylabel('y')
plt.ylim([-4, 4])

plt.xlabel('x')
plt.xlim([0, 80])

colors = np.linspace(0, 1, len(patches))
collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=0.3)
collection.set_array(colors)

ax.add_collection(collection)

plt.show()


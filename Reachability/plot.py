import numpy as np
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

f = open("output2.txt").readlines()

for line in f:
    s = line.split(" ")
    n = int(s[0])
    x = []
    y = []
    for j in range(n):
        px = float(s[2*j+1])
        py = float(s[2*j+2])
        x.append(px)
        y.append(py)

    x.append(x[0])
    y.append(y[0])

    plt.plot(x, y)

plt.show()

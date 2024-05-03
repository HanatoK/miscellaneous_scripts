#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def hilbert_index(x, y, z, level):
    index = 0
    mask = 1 << (level - 1)
    si = 0
    # from https://github.com/AstroUGent/shadowfax/blob/b467282b5d9083a6eaa3a926cf5a638feba08451/src/utilities/Hilbert.cpp#L33C1-L45C76
    hilbert_map = [[( 5,  0), ( 1,  7), ( 4,  1), ( 2,  6), ( 3,  3), ( 3,  4), ( 4,  2), ( 2,  5)],
                   [( 6,  4), ( 2,  7), ( 6,  3), ( 8,  0), ( 0,  5), ( 0,  6), ( 7,  2), ( 7,  1)],
                   [( 1,  6), ( 0,  7), ( 1,  5), ( 9,  4), (10,  1), (11,  0), (10,  2), ( 9,  3)],
                   [( 9,  2), ( 8,  5), ( 0,  3), ( 0,  4), ( 9,  1), ( 8,  6), ( 6,  0), (10,  7)],
                   [( 0,  0), ( 5,  1), ( 8,  3), ( 5,  2), (11,  7), ( 6,  6), ( 8,  4), ( 6,  5)],
                   [( 4,  0), (10,  3), ( 9,  7), (10,  4), ( 0,  1), ( 0,  2), ( 7,  6), ( 7,  5)],
                   [(11,  6), (11,  5), ( 3,  1), ( 3,  2), ( 4,  7), ( 1,  4), ( 9,  0), ( 1,  3)],
                   [( 9,  6), ( 8,  1), ( 5,  7), ( 1,  0), ( 9,  5), ( 8,  2), (11,  4), (11,  3)],
                   [( 1,  2), ( 4,  3), ( 1,  1), ( 7,  0), (10,  5), ( 4,  4), (10,  6), ( 3,  7)],
                   [( 2,  4), ( 5,  5), ( 7,  7), ( 5,  6), ( 2,  3), ( 6,  2), ( 3,  0), ( 6,  1)],
                   [(11,  2), (11,  1), ( 3,  5), ( 3,  6), ( 5,  3), ( 2,  0), ( 5,  4), ( 8,  7)],
                   [( 7,  4), ( 7,  3), ( 4,  5), ( 2,  2), ( 6,  7), (10,  0), ( 4,  6), ( 2,  1)]]
    for i in range(level):
        index <<= 3
        ix = (x & mask) > 0
        iy = (y & mask) > 0
        iz = (z & mask) > 0
        ci = (ix << 2) | (iy << 1) | (iz << 0)
        index |= hilbert_map[si][ci][1]
        si = hilbert_map[si][ci][0]
        mask >>= 1
    return index


level = 3
N = 2**level
cmap = mpl.colormaps['viridis']
point_list = []
for i in range(N):
    for j in range(N):
        for k in range(N):
            point_list.append(np.array([i, j, k, hilbert_index(i, j, k, level)]))
point_list = np.array(point_list)
sorted_list = point_list[point_list[:, -1].argsort()]
# print(sorted_list)
ax = plt.figure().add_subplot(projection='3d')
colors = cmap(sorted_list[:, 3]/len(sorted_list))
x = sorted_list[:, 0]
y = sorted_list[:, 1]
z = sorted_list[:, 2]
points = np.array([x, y, z]).T.reshape(-1, 1, 3)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = Line3DCollection(segments, colors=colors)
ax.add_collection(lc)
ax.set_xlim(x.min(), x.max())
ax.set_ylim(y.min(), y.max())
ax.set_zlim(z.min(), z.max())
plt.show()

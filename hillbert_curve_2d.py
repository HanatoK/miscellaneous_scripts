#!/usr/bin/env python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('QtAgg')

# https://bertvandenbroucke.netlify.app/2019/01/18/space-filling-curves/

def hilbert_index(x, y, level):
    index = 0
    hilbert_map = [[(1, 0), (0, 1), (2, 3), (0, 2)],
                   [(0, 0), (3, 3), (1, 1), (1, 2)],
                   [(2, 2), (2, 1), (0, 3), (3, 0)],
                   [(3, 2), (1, 3), (3, 1), (2, 0)]]
    mask = 1 << (level - 1)
    si = 0
    for i in range(level):
        index <<= 2
        ix = (x & mask) > 0
        iy = (y & mask) > 0
        ci = (ix << 1) | iy
        index |= hilbert_map[si][ci][1]
        si = hilbert_map[si][ci][0]
        mask >>= 1
    return index


# print(hilbert_index(np.uint32(0), np.uint32(1)))

level = 1
N = 2**level
point_list = []
for i in range(N):
    for j in range(N):
        point_list.append(np.array([i, j, hilbert_index(i, j, level)]))
point_list = np.array(point_list)
sorted_list = point_list[point_list[:, -1].argsort()]
print(sorted_list)
plt.plot(sorted_list[:, 0], sorted_list[:, 1])
plt.show()

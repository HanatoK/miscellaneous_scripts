#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

cmap = plt.cm.get_cmap('plasma')

with open('colors.dat', 'w') as fOutput:
    for i in np.linspace(0, 1, 1001):
        for x in cmap(i):
            fOutput.write(f'{x:12.7f}')
        fOutput.write('\n')
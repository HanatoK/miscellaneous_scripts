#!/usr/bin/env python3
import matplotlib
import os
import argparse
matplotlib.use("pgf")
import matplotlib.pyplot as plt
#plt.style.use('fivethirtyeight')
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": False,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 20,
    "axes.linewidth": 2.0,
    "font.size": 16,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{Arimo}",
    ])
})
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

labels = ['Single-Node \n (2 threads)', 'PME GPU \n (24 threads)', 'PME CPU \n (24 threads)', 'PME GPU \n (48 threads)', 'PME CPU \n (48 threads)']
grids_2 = [0.01586297, 0.028019117, 0.0348807, 0.026160833, 0.0255804]
grids_3 = [0.01652817, 0.03444372, 0.04961698, 0.03150375, 0.033056683]
grids_4 = [0.0171544, 0.04079493, 0.04206687, 0.037252133, 0.029302533]
grids_5 = [0.01778213, 0.0474708, 0.05724277, 0.04356487, 0.038752683]
grids_2 = [1.0 / i for i in grids_2]
grids_3 = [1.0 / i for i in grids_3]
grids_4 = [1.0 / i for i in grids_4]
grids_5 = [1.0 / i for i in grids_5]
print(grids_2)
print(grids_3)
print(grids_4)
print(grids_5)

w, h = figaspect(1/1.6)
fig = plt.figure(figsize = (w, h))
width = 0.2
r1 = np.arange(len(labels))
r2 = [i + width for i in r1]
r3 = [i + width for i in r2]
r4 = [i + width for i in r3]
rects1 = plt.barh(r1, grids_2, width, label = '2 Grids')
rects2 = plt.barh(r2, grids_3, width, label = '3 Grids')
rects3 = plt.barh(r3, grids_4, width, label = '4 Grids')
rects4 = plt.barh(r4, grids_5, width, label = '5 Grids')
plt.title('Benchmark on DGX2')
ax = plt.gca()
ax.set_xlabel('Nanoseconds per day')
ax.set_yticks([i + width * 1.5 for i in range(len(labels))])
ax.set_yticklabels(labels)
ax.get_yticklabels()[2].set_color("dimgrey")
ax.get_yticklabels()[4].set_color("dimgrey")
ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False
plt.legend(prop = {'size': 16}, fancybox = False, frameon = False)
plt.grid(axis = 'x', which = 'major', linestyle = '-.')
#plt.grid(axis = 'x', which = 'minor', linestyle = '-.')
ax.set_xlim(0, 76)
y_shift = 0.08
text_size = 12
for i, v in enumerate(grids_2):
    print(i, v)
    percent = (v - grids_2[-1]) / grids_2[-1] * 100
    if percent < 0:
        labelstr = f'{v:.1f}(-{percent:.1f}%)'
    elif percent == 0:
        labelstr = f'{v:.1f}({percent:.1f}%)'
    else:
        labelstr = f'{v:.1f}(+{percent:.1f}%)'
    ax.text(v + 0.05, i - y_shift, labelstr, color='black', fontsize = text_size)
for i, v in enumerate(grids_3):
    print(i, v)
    percent = (v - grids_3[-1]) / grids_3[-1] * 100
    if percent < 0:
        labelstr = f'{v:.1f}(-{percent:.1f}%)'
    elif percent == 0:
        labelstr = f'{v:.1f}({percent:.1f}%)'
    else:
        labelstr = f'{v:.1f}(+{percent:.1f}%)'
    ax.text(v + 0.05, i + width - y_shift, labelstr, color='black', fontsize = text_size)
for i, v in enumerate(grids_4):
    print(i, v)
    percent = (v - grids_4[-1]) / grids_4[-1] * 100
    if percent < 0:
        labelstr = f'{v:.1f}(-{percent:.1f}%)'
    elif percent == 0:
        labelstr = f'{v:.1f}({percent:.1f}%)'
    else:
        labelstr = f'{v:.1f}(+{percent:.1f}%)'
    ax.text(v + 0.05, i + width * 2 - y_shift, labelstr, color='black', fontsize = text_size)
for i, v in enumerate(grids_5):
    print(i, v)
    percent = (v - grids_5[-1]) / grids_5[-1] * 100
    if percent < 0:
        labelstr = f'{v:.1f}(-{percent:.1f}%)'
    elif percent == 0:
        labelstr = f'{v:.1f}({percent:.1f}%)'
    else:
        labelstr = f'{v:.1f}(+{percent:.1f}%)'
    ax.text(v + 0.05, i + width * 3 - y_shift, labelstr, color='black', fontsize = text_size)
fig.tight_layout(pad = 0.2)
plt.savefig('b.png', dpi = 300)

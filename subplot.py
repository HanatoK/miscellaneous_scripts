#!/usr/bin/env python3
import matplotlib
import os
import argparse
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 24,
    "axes.linewidth": 2.0,
    "font.size": 22,
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

# import PMF data
wtm_x, wtm_y = np.genfromtxt('WTM-eABF.dat', unpack = True)
meta_x, meta_y = np.genfromtxt('Meta-eABF.dat', unpack = True)
us_x, us_y = np.genfromtxt('rmsd_US.dat', unpack = True)
us_x = np.arange(45, len(us_y) * 45 + 45, 45) / 1000
wtm_x = np.arange(1, len(wtm_y) + 1, 1) * 3e-4
meta_x = np.arange(1, len(meta_y) + 1, 1) * 3e-4
X = [us_x, wtm_x, meta_x]
Y = [us_y, wtm_y, meta_y]
# setup labels and colors
labels = ['US/WHAM', 'WTM-eABF', 'Meta-eABF\n (small Gaussian)']
colors = ['green', 'royalblue', 'red']
alpha_ratio = 0.6

w, h = figaspect(1/2.2)
fig = plt.figure(figsize = (w,h))
# first subplot
plt.subplot(1,2,1)
plt.plot(X[0], Y[0], label = labels[0], color = colors[0], alpha = alpha_ratio)
plt.plot(X[1], Y[1], label = labels[1], color = colors[1], alpha = alpha_ratio)
plt.plot(X[2], Y[2], label = labels[2], color = colors[2], alpha = alpha_ratio)
ax = plt.gca()
ax.set_xlim(0 - 0.02, 0.2+0.05)
ax.set_ylim(-1, 50)
ax.xaxis.set_major_locator(MultipleLocator(0.05))
ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False

# second subplot
plt.subplot(1,2,2)
plt.plot(X[0], Y[0], label = labels[0], color = colors[0], alpha = alpha_ratio)
plt.plot(X[1], Y[1], label = labels[1], color = colors[1], alpha = alpha_ratio)
plt.plot(X[2], Y[2], label = labels[2], color = colors[2], alpha = alpha_ratio)
ax = plt.gca()
num_ticks = 6
ax.set_xlim(0.25, 1.5)
ax.set_ylim(-0.05, 1.25)
ax.xaxis.set_major_locator(plt.MaxNLocator(num_ticks))
ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False

#plt.subplots_adjust(wspace = 0.25)
plt.tight_layout(pad = 0.1)
plt.savefig('rmsd_split.png', dpi=600, transparent=False)

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
from matplotlib.ticker import AutoMinorLocator

def plot_rmsd(dataset_x, dataset_y, labels, colors, outputname):
    w, h = figaspect(1/1.1)
    plt.figure(figsize = (w,h))
    alpha_ratio = 1.0
    for x, y, label, color in zip(dataset_x, dataset_y, labels, colors):
        plt.plot(x, y, label = label, color = color, alpha = alpha_ratio, linewidth = 1.0)
    ax = plt.gca()
    plt.xlabel(r'Time (ns)')
    plt.ylabel('RMSD($\Delta G$) (kcal/mol)')
    ax.set_xlim(0, 100)
    ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
    ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.legend(prop = {'size': 14}, fancybox = False, frameon = False)
    plt.tight_layout(pad = 0.1)
    plt.savefig(outputname, dpi=600, transparent=False)
    plt.close()

x1, y1 = np.genfromtxt('rmsd_1000K.dat', unpack = True)
x2, y2 = np.genfromtxt('rmsd_2000K.dat', unpack = True)
x3, y3 = np.genfromtxt('rmsd_4000K.dat', unpack = True)
x4, y4 = np.genfromtxt('rmsd_8000K.dat', unpack = True)
x5, y5 = np.genfromtxt('rmsd_infinite.dat', unpack = True)
labels = ['1000 K', '2000 K', '4000 K', '8000 K', 'Infinite']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
time_ns = np.linspace(1, len(y1), len(y1))
dataset_x = [time_ns, time_ns, time_ns, time_ns, time_ns]
dataset_y = [y1, y2, y3, y4, y5]
plot_rmsd(dataset_x, dataset_y, labels, colors, 'rmsd_temperature.png')

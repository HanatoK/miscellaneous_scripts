#!/usr/bin/python3
import argparse
import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "Arimo",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 16,
    "axes.linewidth": 1.5,
    "font.size": 12,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",  # serif font via preamble
    ])
})
import numpy as np
from matplotlib import ticker
import math
from matplotlib.figure import figaspect
from matplotlib.patches import Rectangle

x1, y1, z1 = np.genfromtxt('hist_elect.dat', unpack=True)
x2, y2, z2 = np.genfromtxt('hist_vdw.dat', unpack=True)
x3, y3, z3 = np.genfromtxt('hist_intra_elect.dat', unpack=True)
x4, y4, z4 = np.genfromtxt('hist_intra_vdw.dat', unpack=True)

X = [x1, x2, x3, x4]
Y = [y1, y2, y3, y4]
Z = [z1, z2, z3, z4]
titles = ['Electrostatic', 'VdW', 'Electrostatic', 'VdW']
fig_labels = ['A', 'B', 'C', 'D']

def reshape_xyz(x, y, z):
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    return xi, yi, zi


w, h = figaspect(1/1.2)
fig, axs = plt.subplots(2, 2, figsize=(w, h), sharex=True, sharey=True)

for row in range(2):
    for col in range(2):
        i = row * 2 + col
        xi, yi, zi = reshape_xyz(X[i], Y[i], Z[i])
        cf = axs[row, col].contourf(xi, yi, zi, levels=40, cmap=plt.cm.turbo, zorder=2)
        clb = fig.colorbar(cf, ax=axs[row, col])
        axs[row, col].set_title(titles[i], fontsize=16)
        axs[row, col].text(3, 1, 'Basin 1', fontsize=10, zorder=5)
        axs[row, col].text(-2.5, 8, 'Basin 2', fontsize=10, zorder=5)
        axs[row, col].text(5, 14.25, 'Basin 3', fontsize=10, zorder=5)
        axs[row, col].text(1, 16, fig_labels[i], fontsize=28, fontweight='bold', zorder=5)
        clb.ax.tick_params(labelsize=12)
        clb.ax.set_title("kcal/mol", fontsize=10)
        clbticks = [float(i.get_position()[1]) for i in clb.ax.get_yticklabels()]
        clbticksstr = ['{:.0f}'.format(i) for i in clbticks]
        print(clbticksstr)
        clb.ax.set_yticklabels(clbticksstr, fontsize = 10)


for row in range(2):
    for col in range(2):
        axs[row, col].tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True, labelsize=16)
        axs[row, col].tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True, labelsize=16)
        axs[row, col].xaxis.get_major_formatter()._usetex = False
        axs[row, col].yaxis.get_major_formatter()._usetex = False
        axs[row, col].grid(zorder=1)
        axs[row, col].add_patch(Rectangle((3, 2.5), width = 5, height = 1.5, fill = None, linewidth = 1.0, color = 'black', zorder=3))
        axs[row, col].add_patch(Rectangle((2.5, 5), width = 1.5, height = 3, fill = None, linewidth = 1.0, color = 'black', zorder=3))
        axs[row, col].add_patch(Rectangle((5, 8), width = 7, height = 6, fill = None, linewidth = 1.0, color = 'black', zorder=3))
        axs[row, col].set_xlim(0, 20)
        axs[row, col].set_ylim(0, 20)
        axs[row, col].xaxis.set_major_locator(plt.MultipleLocator(5))
        axs[row, col].yaxis.set_major_locator(plt.MultipleLocator(5))

#fig.add_subplot(111, frameon=False)
#plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
#plt.grid(False)
#plt.xlabel(r'$d_1$ (Asp3N-Gly7O)')
#plt.ylabel(r'$d_2$ (Asp3N-Thr8O)')
fig.text(0.5, 0.01, r'$d_1$ (Asp3N-Gly7O) (Å)', ha='center', fontsize=20)
fig.text(0.01, 0.5, r'$d_2$ (Asp3N-Thr8O) (Å)', va='center', rotation='vertical', fontsize=20)
#plt.tight_layout(pad = 0.1)
plt.savefig('all_pair_interaction.png', dpi=400, bbox_inches = 'tight', transparent=False)
plt.savefig('all_pair_interaction.pdf', dpi=400, bbox_inches = 'tight', transparent=False)

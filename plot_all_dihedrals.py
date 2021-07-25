#!/usr/bin/env python3
import pandas as pd
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",
    ])
})
from matplotlib.figure import figaspect
from matplotlib.ticker import AutoMinorLocator

#import data
dihedrals = pd.read_csv('path-reparam-me.txt', delimiter='\s+', header=None)

w, h = figaspect(1/1.5)
plt.figure(figsize = (w,h))
X = list(range(1, len(dihedrals)+1))
plt.plot(X, dihedrals[0], label='$\phi_1$', color='orangered')
plt.plot(X, dihedrals[1], label='$\phi_2$', color='brown')
plt.plot(X, dihedrals[2], label='$\phi_3$', color='darkorange')
plt.plot(X, dihedrals[3], label='$\psi_1$', color='limegreen')
plt.plot(X, dihedrals[4], label='$\psi_2$', color='darkcyan')
plt.plot(X, dihedrals[5], label='$\psi_3$', color='royalblue')
plt.scatter(X, dihedrals[0], marker='*', color='orangered')
plt.scatter(X, dihedrals[1], marker='*', color='brown')
plt.scatter(X, dihedrals[2], marker='*', color='darkorange')
plt.scatter(X, dihedrals[3], marker='*', color='limegreen')
plt.scatter(X, dihedrals[4], marker='*', color='darkcyan')
plt.scatter(X, dihedrals[5], marker='*', color='royalblue')
plt.xlabel('Image index')
plt.ylabel('Dihedral angle (°)')
ax = plt.gca()
ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False
ax.set_xlim(0, 35)
ax.set_ylim(-180, 180)
ax.xaxis.set_major_locator(plt.MultipleLocator(5))
ax.yaxis.set_major_locator(plt.MultipleLocator(60))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.legend(prop = {'size': 22}, fancybox = False, frameon = False, ncol=3)
plt.savefig('all.png', dpi=300, bbox_inches='tight')
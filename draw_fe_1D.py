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
    "axes.labelsize": 20,
    "axes.linewidth": 2.0,
    "font.size": 16,
    "pgf.preamble": [
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{Arimo}",  # serif font via preamble
         ]
})
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

parser = argparse.ArgumentParser()
parser.add_argument("pmf", help = "specify the PMF file")
parser.add_argument("-o", "--output", help = "specify the PNG output image file")
parser.add_argument("-c", "--colvars", action = "store_true", help = "use Colvars units")
parser.add_argument("-s", "--save", action = "store_true", help = "save new PMF file with min-to-zero")
args = parser.parse_args()
pmf = args.pmf
if args.output is None:
    png = pmf + ".png"
else:
    png = args.output

x, y = np.genfromtxt(pmf, unpack=True)
y = y - min(y)
w, h = figaspect(1/1.1)
plt.figure(figsize = (w,h))
if args.colvars == False:
    y = y / 4.184
plt.plot(x, y)
#if args.colvars == True:
    #plt.xlabel("Distance(Å)")
#else:
    #plt.xlabel("Distance(nm)")
plt.xlabel(r"$\xi (s)$")
plt.ylabel(r"$\Delta G$ (kcal/mol)")
ax = plt.gca()
ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False
ax.xaxis.set_major_locator(plt.MaxNLocator(10, prune = 'lower'))
ax.yaxis.set_major_locator(plt.MaxNLocator(10))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.savefig(png, dpi=400, bbox_inches = 'tight', transparent=True)

if args.save is True:
    fn = os.path.basename(pmf)
    newpmf = "new_" + fn
    np.savetxt(newpmf, np.transpose([x, y]), delimiter = "\t", fmt = "%10.7f")

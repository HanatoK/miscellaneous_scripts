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
    "font.size": 20,
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

parser = argparse.ArgumentParser()
parser.add_argument("grad", help = "specify the gradient file")
parser.add_argument("-o", "--output", help = "specify the output prefix")
args = parser.parse_args()

def integrate1D(x, grad):
    width = x[1] - x[0]
    start_x = x[0] - 0.5 * width
    start_pmf = 0
    pmf_cv = [start_x]
    pmf = [start_pmf]
    for g_i in np.nditer(grad):
        start_x += width
        start_pmf += g_i * width
        pmf_cv.append(start_x)
        pmf.append(start_pmf)
    return np.array(pmf_cv), np.array(pmf)

def plot_fes_1D(cv, pmf, xtitle, ytitle, label, color, outputname):
    w, h = figaspect(1/1.1)
    plt.figure(figsize = (w,h))
    alpha_ratio = 1.0
    plt.plot(cv, pmf, label = label, color = color, alpha = alpha_ratio)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    ax = plt.gca()
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
    ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
    plt.legend(prop = {'size': 14}, fancybox = False, frameon = False)
    plt.tight_layout()
    plt.savefig(outputname, dpi=600, transparent=False)

x, grad_x = np.genfromtxt(args.grad, unpack = True)
cv, pmf = integrate1D(x, grad_x)
pmf = pmf - np.min(pmf)
outputprefix = args.output
outputpmf = outputprefix + '.pmf'
outputpng = outputprefix + '.png'
np.savetxt(outputpmf, np.c_[cv, pmf], fmt = '%10.2f %12.7f')
plot_fes_1D(cv, pmf, r'$i$', '$\Delta G$ (kcal/mol)', None, 'red', outputpng)

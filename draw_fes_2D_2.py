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
parser.add_argument("-o", "--output", default = "output.png", help = "specify the PNG output image file")
parser.add_argument("--xtitle", default = "CV1", help = "title along X axis")
parser.add_argument("--ytitle", default = "CV2", help = "title along Y axis")
parser.add_argument("--levels", default = 25, type = int, help = "number of levels")
args = parser.parse_args()

def plotfes(pmffilename, pngfilename, xtitle, ytitle):
    x, y, z = np.genfromtxt(pmffilename, unpack=True)
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    fig = plt.figure()
    cf = plt.contourf(xi, yi, zi, args.levels, cmap=plt.cm.jet)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title('Free Energy Surface')
    ax = plt.gca()
    ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
    ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    #ax.xaxis.set_major_locator(MultipleLocator(0.01))
    #ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    clb = plt.colorbar()
    clb.ax.set_title("kcal/mol")
    clb.ax.xaxis.get_major_formatter()._usetex = False
    clb.ax.yaxis.get_major_formatter()._usetex = False
    clbticks = list(cf.levels)
    clbticksstr = ['{:.1f}'.format(i) for i in clbticks]
    clb.ax.set_title("Error (kcal/mol)", fontsize = 16, pad = 10.0)
    clb.ax.set_yticklabels(clbticksstr, fontsize = 16)
    plt.tight_layout(pad = 0.1)
    plt.savefig(pngfilename, dpi=400, transparent=True)
    return

plotfes(args.pmf, args.output, xtitle = args.xtitle, ytitle = args.ytitle)

#!/usr/bin/env python3
import matplotlib
import math
import argparse
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
    "axes.unicode_minus": False,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",
    ])
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
parser.add_argument("-t", "--title", default="Free Energy Surface", help="title of figure")
args = parser.parse_args()
pmf = args.pmf
colvars = False
if args.colvars:
    colvars = True
if args.output is None:
    png = pmf + ".png"
else:
    png = args.output

def plotfes(pmffilename, pngfilename, title='Free Energy Surface'):
    w, h = figaspect(1/1)
    plt.figure(figsize = (w,h))
    #px1, py1 = np.genfromtxt("path_C36.txt", unpack=True)
    #px2, py2 = np.genfromtxt("path.txt", unpack=True)
    x, y, z = np.genfromtxt(pmffilename, unpack=True)
    # z = np.clip(z, np.min(z), 20)
    if colvars is False:
        x = x / math.pi * 180
        y = y / math.pi * 180
        z = z / 4.184
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    z_range = np.linspace(0, 7.2, 25)
    plt.figure()
    cf = plt.contourf(xi, yi, zi, levels=z_range, cmap='turbo')
    #plt.plot(px1, py1, color = 'orange')
    #plt.scatter(px1, py1, color = 'grey')
    #plt.plot(px2, py2, color = 'red')
    #plt.scatter(px2, py2, color = 'black')
    ax = plt.gca()
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.tick_params(direction='in', which='major', length=6.0,
                   width=1.0, top=True, right=True, pad=8.0)
    ax.tick_params(direction='in', which='minor', length=3.0,
                   width=1.0, top=True, right=True, pad=8.0)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.xaxis.set_major_locator(plt.MultipleLocator(60))
    ax.yaxis.set_major_locator(plt.MultipleLocator(60))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel(r'$\phi$ (degree)')
    plt.ylabel(r'$\psi$ (degree)')
    plt.title(title)
    clb = plt.colorbar()

    clb.ax.set_title("kcal/mol", pad=10.0, fontsize=22)
    clb.ax.xaxis.get_major_formatter()._usetex = False
    clb.ax.yaxis.get_major_formatter()._usetex = False
    # clbticks = [i.get_text() for i in clb.ax.get_yticklabels()]
    # clbticksstr = [i.strip('$') for i in clbticks]
    # print(clbticks)
    # print(clbticksstr)
    # clb.ax.set_yticklabels(clbticksstr, fontsize = 22)
    clb.ax.tick_params(labelsize=22)
    # plt.tight_layout(pad = 0.2)
    plt.savefig(pngfilename, dpi=400, bbox_inches='tight', transparent=False)
    return

plotfes(pmf, png, args.title)

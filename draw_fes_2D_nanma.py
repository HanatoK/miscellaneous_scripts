#!/usr/bin/python3
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({
    #"pgf.texsystem": "lualatex",
    "font.family": "Arimo",  # use serif/main font for text elements
    #"text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 16,
    "axes.linewidth": 2.0,
    "font.size": 12,
    #"pgf.preamble": [
         #"\\usepackage{units}",          # load additional packages
         #"\\usepackage{metalogo}",
         #"\\usepackage{unicode-math}",   # unicode math setup
         #r"\setmathfont{MathJax_Math}",
         #r"\setmainfont{FreeSans}",  # serif font via preamble
         #]
})
import numpy as np
import math

parser = argparse.ArgumentParser()
parser.add_argument("pmf", help = "specify the PMF file")
parser.add_argument("-o", "--output", help = "specify the PNG output image file")
parser.add_argument("-c", "--colvars", action = "store_true", help = "use Colvars units")
args = parser.parse_args()
pmf = args.pmf
colvars = False
if args.colvars:
    colvars = True
if args.output is None:
    png = pmf + ".png"
else:
    png = args.output

def plotfes(pmffilename, pngfilename):
    #px1, py1 = np.genfromtxt("path_C36.txt", unpack=True)
    #px2, py2 = np.genfromtxt("path.txt", unpack=True)
    x, y, z = np.genfromtxt(pmffilename, unpack=True)
    z = np.clip(z, np.min(z), 20)
    if colvars is False:
        x = x / math.pi * 180
        y = y / math.pi * 180
        z = z / 4.184
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    plt.figure()
    cf = plt.contourf(xi, yi, zi, 40, cmap=plt.cm.jet)
    #plt.plot(px1, py1, color = 'orange')
    #plt.scatter(px1, py1, color = 'grey')
    #plt.plot(px2, py2, color = 'red')
    #plt.scatter(px2, py2, color = 'black')
    ax = plt.gca()
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.xaxis.set_major_locator(plt.MaxNLocator(10))
    #ax.xaxis.set_minor_locator(plt.MaxNLocator(50))
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    #ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
    plt.xlabel(r'$\phi$ (degree)')
    plt.ylabel(r'$\psi$ (degree)')
    #plt.title('Free Energy Surface')
    clb = plt.colorbar(cf)
    clb.ax.set_title("kcal/mol")
    ticklabs = clb.ax.get_yticklabels()
    clb.ax.set_yticklabels(ticklabs,ha='right')
    clb.ax.yaxis.set_tick_params(pad=25)
    plt.tight_layout()
    plt.savefig(pngfilename, dpi=400, transparent=False)
    return

plotfes(pmf, png)

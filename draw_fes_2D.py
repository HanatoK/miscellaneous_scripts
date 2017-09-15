#!/usr/bin/python3
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    x, y, z = np.genfromtxt(pmffilename, unpack=True)
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
    plt.contourf(xi, yi, zi, 25, cmap=plt.cm.jet)
    plt.xlabel(r'$\phi$(degree)')
    plt.ylabel(r'$\psi$(degree)')
    plt.title('Free Energy Surface')
    clb = plt.colorbar()
    clb.ax.set_title("kCal/mol")
    plt.savefig(pngfilename, dpi=400, transparent=True)
    return

plotfes(pmf, png)

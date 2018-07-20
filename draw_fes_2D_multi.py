#!/usr/bin/env python3
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math

parser = argparse.ArgumentParser()
parser.add_argument('pmf', nargs = '+', help = 'the PMF file')
parser.add_argument('-o', '--output', default = 'output.png', help = 'the PNG output image file')
parser.add_argument('--max', type = float, help = 'upperbound of the FES')
parser.add_argument('--unit', type = float, default = 1.0, help = 'energy unit')
parser.add_argument('--xlabel', default = 'CV1', help = 'label of the x axis')
parser.add_argument('--ylabel', default = 'CV2', help = 'label of the y axis')
parser.add_argument('--fontsize', type = float, default = 9.0, help = 'size of the font')
parser.add_argument('--plumed', action = 'store_true', help = 'is this a PLUMED fes from sum_hills?')
args = parser.parse_args()

'''function to plot PMF'''
def plotPMF(x, y, z, maxpmf, xlabel, ylabel, unit, outfile):
    z = z / unit
    fig = plt.figure()
    colorbarBins = (maxpmf - 0) / 2.0 + 1
    cf = plt.contourf(x, y, z, levels = np.linspace(0, maxpmf, colorbarBins), cmap=plt.cm.jet)
    plt.title('Free Energy Surface')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    clb = plt.colorbar(cf)
    clb.ax.set_title("kcal/mol")
    plt.savefig(outfile, dpi = 400)
    plt.close(fig)

'''update the font size'''
matplotlib.rcParams.update({'font.size': args.fontsize})

for pmfData in args.pmf:
    if args.plumed is True:
        x, y, z, dx, dy = np.genfromtxt(pmfData, unpack = True)
        binx = len(set(x))
        biny = len(set(y))
        xi = x.reshape(biny, binx)
        yi = y.reshape(biny, binx)
        zi = z.reshape(biny, binx)
    else:
        x, y, z = np.genfromtxt(pmfData, unpack = True)
        binx = len(set(x))
        biny = len(set(y))
        xi = x.reshape(binx, biny)
        yi = y.reshape(binx, biny)
        zi = z.reshape(binx, biny)
    if len(args.pmf) == 1:
        pngfile = args.output
    else:
        pngfile = pmfData + '.png'
    if args.max is None:
        maxpmf = max(z)
    else:
        maxpmf = args.max
    plotPMF(xi, yi, zi, maxpmf, args.xlabel, args.ylabel, args.unit, pngfile)

#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib
import os
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("pmf", nargs = '+', help = "specify the PMF file")
parser.add_argument("-o", "--output", help = "specify the PNG output image file")
parser.add_argument("-c", "--colvars", action = "store_true", help = "use Colvars units")
parser.add_argument("-r", "--reference", help = "reference PMF file")
args = parser.parse_args()
png = args.output

if args.reference is not None:
    rx, ry = np.genfromtxt(args.reference, unpack=True)
    plt.plot(rx, ry, label = 'Reference')

for pmffile in args.pmf:
    x, y = np.genfromtxt(pmffile, unpack=True)
    y = y + abs(min(y))
    if args.colvars == False:
        y = y / 4.184
    plt.plot(x, y, label = pmffile)
plt.ylabel("Free energy(kcal/mol)")
ax = plt.gca()
ax.tick_params(direction="in")
plt.legend(loc = 'upper left', fontsize=7, labelspacing=0.5)
plt.savefig(png, dpi=400, transparent=True)

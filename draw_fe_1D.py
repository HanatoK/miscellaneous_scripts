#!/usr/bin/python3
import argparse
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("pmf", help = "specify the PMF file")
parser.add_argument("-o", "--output", help = "specify the PNG output image file")
parser.add_argument("-c", "--colvars", action = "store_true", help = "use Colvars units")
args = parser.parse_args()
pmf = args.pmf
if args.output is None:
    png = pmf + ".png"
else:
    png = args.output

x, y = np.genfromtxt(pmf, unpack=True)
y = y + abs(min(y))
if args.colvars == False:
    y = y / 4.184
plt.plot(x, y)
plt.xlabel("Distance(nm)")
plt.ylabel("Free energy(kCal/mol)")
ax = plt.gca()
ax.tick_params(direction="in")
plt.savefig(png, dpi=400, transparent=True)

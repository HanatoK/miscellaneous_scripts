#!/usr/bin/python3
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('rmsd', nargs = '+', help = 'specify the RMSD data file(s)')
parser.add_argument('-o', '--output', help = 'PNG output file')
parser.add_argument('-l', '--line', type = int, help = 'PNG output file')
args = parser.parse_args()
data = args.rmsd
outpng = args.output


plt.figure()
y = np.genfromtxt(data[0], unpack=True, usecols=1)
if args.line is None:
    nl = len(y)
else:
    nl = args.line
x = np.linspace(1, nl, nl)
for file in data:
    lb = os.path.splitext(file)[0]
    y = np.genfromtxt(file, unpack=True, usecols=1)
    y = y[0:nl]
    plt.plot(x,y,label=lb)
plt.xlabel('Simulation time(ns)')
plt.ylabel('Free energy RMSD(kcal/mol)')
plt.tick_params(direction='in')
leg = plt.legend()
leg.get_frame().set_alpha(0.0)
plt.savefig(outpng, dpi=400, transparent=True)

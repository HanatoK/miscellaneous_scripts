#!/usr/bin/env python3
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import gridspec

parser = argparse.ArgumentParser()
parser.add_argument("traj", help = "specify the collective variables trajectory file")
parser.add_argument("-c", "--column", type = int, default = 1, help = "specify the column of CV in traj file")
parser.add_argument("-o", "--output", default = "output", help = "specify the output image")
parser.add_argument("-n", "--frequency", type = int, default = 1000, help = "how many CV data points in a plot")
args = parser.parse_args()

traj = args.traj
col = args.column
outfile = args.output
freq = args.frequency

data = np.genfromtxt(traj, unpack = True, )
cvlength = len(data[col-1])
maxcv = max(data[col-1])
mincv = min(data[col-1])
step = data[0]

def plotCV(start, end, cv, filename):
    y = cv[start:end]
    fig = plt.figure()
    plt.plot(step[start:end], y, linewidth = 0.5, color = "red")
    plt.ylabel("CV")
    axes = plt.gca()
    axes.set_ylim([mincv, maxcv])
    axes.tick_params(direction="in")
    axes.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    outpng = filename + ".png"
    plt.savefig(outpng, dpi=400, transparent=False)
    plt.close(fig)

def subplotCV(cv, filename):
    plotColumns = 3
    nplots = int(cvlength / freq)
    plotRows = int(math.ceil(nplots / plotColumns))
    gs = gridspec.GridSpec(plotRows, plotColumns)
    fig = plt.figure()
    for i in range(1, nplots + 1):
        ax = fig.add_subplot(gs[i-1])
        ax.tick_params(direction="in")
        ax.tick_params(axis='x', labelbottom=False)
        #ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.set_ylim([mincv, maxcv])
        start = freq * (i - 1)
        end = freq * i
        y = cv[start:end]
        ax.plot(step[start:end], y, linewidth = 0.25, color = "red")
    ax = fig.add_subplot(gs[nplots])
    ax.tick_params(direction="in")
    ax.tick_params(axis='x', labelbottom=False)
    #ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.set_ylim([mincv, maxcv])
    start = freq * nplots
    end = cvlength
    y = cv[start:end]
    ax.plot(step[start:end], y, linewidth = 0.25, color = "red")
    outpng = filename + ".png"
    plt.savefig(outpng, dpi=400, transparent=False)
    plt.close(fig)
        

loop_count = int(cvlength / freq)
for i in range(1, loop_count + 1):
    start = freq * (i - 1)
    end = freq * i
    filename = outfile + str(i)
    plotCV(start, end, data[col-1], filename)

subplotCV(data[col-1], outfile + "_multi")

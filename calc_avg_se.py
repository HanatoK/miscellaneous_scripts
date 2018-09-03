#!/usr/bin/env python3
import numpy as np
import matplotlib
import os
matplotlib.use("Agg")
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

def calc_avg_se(datalist):
    sum1 = None
    sum2 = None
    i = 0
    for l in datalist:
        if (i == 0):
            sum1 = l
            sum2 = l * l
        else:
            sum1 += l
            sum2 += l * l
        i = i + 1
    avg1 = sum1 / i
    avg2 = sum2 / i
    se = np.sqrt(avg2 - avg1 * avg1) / np.sqrt(i)
    sd = np.sqrt(avg2 - avg1 * avg1)
    return [avg1, se, sd]
    

x1, y1 = np.genfromtxt("output_best.pmf", unpack=True)
x2, y2 = np.genfromtxt("output_3.pmf", unpack=True)
x3, y3 = np.genfromtxt("output_4.pmf", unpack=True)
ly = [y1, y2, y3]
avg, se, sd = calc_avg_se(ly)
np.savetxt("stat.pmf", np.c_[x1, avg, sd], fmt="%.7f")

plt.figure()
plt.ylabel("Free energy (kcal/mol)")
plt.xlabel("End-to-end distance (Å)")
ax = plt.gca()
ax.tick_params(direction="in")
plt.errorbar(x1, avg, yerr=sd, ecolor="lightcoral", color="red", label="eABF+GaMD(only dihedral boost)")
rx, ry = np.genfromtxt("ref.dat", unpack=True)
plt.plot(rx, ry, label="Reference")
plt.legend(loc = 'upper left', fontsize=11, labelspacing=0.5, fancybox=True, framealpha=0)
plt.savefig("stat.png", dpi=400, transparent=True)

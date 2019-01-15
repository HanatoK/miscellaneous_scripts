#!/usr/bin/env python3
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import copy

'''boltzmann constant in kcal/mol'''
kb = 0.0019872041
'''temperature in K'''
t = 300
parser = argparse.ArgumentParser()
parser.add_argument("histo", help = "specify the histogram file")
parser.add_argument("-o", "--output", help = "specify output prefix")
args = parser.parse_args()
histo = args.histo
outpmf = args.output + ".pmf"
outpng = args.output + ".png"

'''function for calculating pmf'''
def calcpmf(probability):
    '''copy'''
    tmp = copy.deepcopy(probability)
    '''allow modify np array'''
    total = np.sum(probability)
    for p_i in np.nditer(tmp, op_flags=['readwrite']):
        if p_i > 0:
            p_i[...] = - kb * t * np.log(p_i / total)
    Amax = max(tmp)
    print("maximum PMF:", Amax)
    for p_i in np.nditer(probability, op_flags=['readwrite']):
        if p_i > 0:
            p_i[...] = - kb * t * np.log(p_i / total)
        else:
            p_i[...] = Amax
    Amin = min(probability)
    for p_i in np.nditer(probability, op_flags=['readwrite']):
        p_i[...] -= Amin
    return probability

def plotpmf(pngfilename,x,z):
    plt.plot(x, z)
    ax = plt.gca()
    ax.tick_params(direction="in")
    plt.savefig(pngfilename, dpi=400, transparent=True)

x, z = np.genfromtxt(histo, unpack=True)
z = calcpmf(z)
fpmf = open(outpmf, "w")
for i, j in zip(x, z):
    print("%.7f   %.7f"%(i,j),file=fpmf)
fpmf.close()
plotpmf(outpng,x,z)

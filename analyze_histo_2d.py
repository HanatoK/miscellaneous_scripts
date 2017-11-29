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
    for p_i in np.nditer(tmp, op_flags=['readwrite']):
        if p_i > 0:
            p_i[...] = - kb * t * np.log(p_i)
    Amax = max(tmp)
    print("maximum PMF:", Amax)
    for p_i in np.nditer(probability, op_flags=['readwrite']):
        if p_i > 0:
            p_i[...] = - kb * t * np.log(p_i)
        else:
            p_i[...] = Amax
    return probability

x, y, z = np.genfromtxt(histo, unpack=True)
z = calcpmf(z)

'''function for 2D pmf plotting'''
def plotpmf(pngfilename,x,y,z):
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
    clb.ax.set_title("kcal/mol")
    plt.savefig(pngfilename, dpi=400, transparent=True)
    return

'''write PMF to file'''
fpmf = open(outpmf, "w")
for i,j,k in zip(x,y,z):
    print("%.7f   %.7f   %.7f"%(i,j,k),file=fpmf)
fpmf.close()
plotpmf(outpng,x,y,z)

#!/usr/bin/env python3
#import matplotlib
#import os
import argparse
#matplotlib.use("pgf")
#import matplotlib.pyplot as plt
#plt.rcParams.update({
    #"pgf.texsystem": "lualatex",
    #"font.family": "serif",  # use serif/main font for text elements
    #"text.usetex": True,     # use inline math for ticks
    #"pgf.rcfonts": False,    # don't setup fonts from rc parameters
    #"axes.labelsize": 24,
    #"axes.linewidth": 2.0,
    #"font.size": 22,
    #"pgf.preamble": '\n'.join([
         #"\\usepackage{units}",
         #"\\usepackage{metalogo}",
         #"\\usepackage{unicode-math}",
         #r"\setmathfont{MathJax_Math}",
         #r"\setmainfont{Arimo}",
    #])
#})
import numpy as np
#from matplotlib.figure import figaspect
#from scipy.interpolate import interp1d
#import matplotlib.ticker as ticker
#from matplotlib.ticker import AutoMinorLocator

def compute_derivatives_1D(cv, pmf):
    width = cv[1] - cv[0]
    x = [cv[0] + width * 0.5]
    grad = [(pmf[1] - pmf[0]) / width]
    for i in range(1, len(cv) - 1):
        x.append(x[-1] + width)
        grad.append((pmf[i + 1] - pmf[i]) / width)
    return x, grad

def generate_dummy_count(cv):
    width = cv[1] - cv[0]
    x = [cv[0] + width * 0.5]
    count = [1]
    for i in range(1, len(cv) - 1):
        x.append(x[-1] + width)
        count.append(1)
    return x, count

def generate_header(cv):
    line1 = '# 1'
    line2 = '# ' + format(cv[0], '.4f') + ' ' + format(cv[1] - cv[0], '.4f') + ' ' + str(len(cv) - 1) + ' 0\n'
    result = line1 + '\n' + line2
    return result

#parser = argparse.ArgumentParser()
#parser.add_argument("pmf", help = "specify the PMF file")
#parser.add_argument("-o", "--output", help = "specify the output prefix")
#args = parser.parse_args()

#cv, pmf = np.genfromtxt(args.pmf, unpack = True)
#x, grad = compute_derivatives_1D(cv, pmf)
#outputgrad = args.output + '.grad'
#np.savetxt(outputgrad, np.c_[x, grad], fmt = '%8.2f %20.15f')
#x, count = generate_dummy_count(cv)
#outputcount = args.output + '.count'
#np.savetxt(outputcount, np.c_[x, count], fmt = '%8.2f %d')

#window_prefix = 'window-1'

for j in range(1, 10):
    window_prefix = 'window-' + str(j)
    for i in range(0, 9009):
        pmf_file = window_prefix + '.hist_' + str(i).zfill(4) + '.pmf'
        grad_name = window_prefix + '.hist_' + str(i).zfill(4) + '.grad'
        count_name = window_prefix + '.hist_' + str(i).zfill(4) + '.count'
        f_grad = open(grad_name, 'w')
        cv, pmf = np.genfromtxt(pmf_file, unpack = True)
        header = generate_header(cv)
        f_grad.write(header)
        x, grad = compute_derivatives_1D(cv, pmf)
        np.savetxt(f_grad, np.c_[x, grad], fmt = '%8.2f %20.15f')
        f_count = open(count_name, 'w')
        f_count.write(header)
        x, count = generate_dummy_count(cv)
        np.savetxt(f_count, np.c_[x, count], fmt = '%8.2f %d')
        f_grad.close()
        f_count.close()

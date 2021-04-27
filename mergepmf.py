#!/usr/bin/env python3
import matplotlib
import os
import argparse
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": False,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 24,
    "axes.linewidth": 2.0,
    "font.size": 22,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{Arimo}",
    ])
})
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import copy

'''boltzmann constant in kcal/mol'''
kb = 0.0019872041
'''temperature in K'''
t = 300

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

def readHeaderString(inputfile):
    headerLines = []
    with open(inputfile, 'r') as fp_in:
        for line in fp_in:
            if line.strip():
                fields = line.split()
                if fields[0][0] == '#':
                    headerLines.append(line)
    return headerLines

def mergepmf(pmf1, pmf2, outputname, plot = True):
    # merge PMFs
    cv_czar, pmf_czar = np.genfromtxt(pmf1, unpack = True)
    cv_amd, pmf_amd = np.genfromtxt(pmf2, unpack = True)
    #cv_count, abf_count = np.genfromtxt('deca.count', unpack = True)
    #abf_count_correction = calcpmf(abf_count)
    #pmf_amd += abf_count_correction
    f = interp1d(cv_amd, pmf_amd, fill_value = 'extrapolate', kind = 'quadratic')
    pmf_amd = f(cv_czar)
    pmf_total = pmf_amd + pmf_czar
    pmf_total_min = np.min(pmf_total)
    pmf_total = pmf_total - pmf_total_min
    output_pmf = outputname + '.pmf'
    headerLines = readHeaderString(pmf1)
    with open(output_pmf, 'w') as fp_out:
        for line in headerLines:
            fp_out.write(line)
        for cv, pmf in np.nditer([cv_czar, pmf_total]):
            fp_out.write(f'{cv:10.4f} {pmf:12.7f}\n')
    #np.savetxt(output_pmf, np.c_[cv_czar, pmf_total], fmt = '%10.4f %12.7f')

    if plot is True:
        # reference PMF
        ref_cv, ref_pmf = np.genfromtxt('../ref.dat', unpack = True)
        # WTM-eABF PMF
        wtm_cv, wtm_pmf = np.genfromtxt('wtm-eabf.pmf', unpack = True)

        # plotting
        w, h = figaspect(1/1.2)
        plt.figure(figsize = (w,h))
        plt.plot(ref_cv, ref_pmf, color = 'black', label = 'Reference')
        plt.plot(wtm_cv, wtm_pmf, color = 'tab:orange', label = 'WTM-eABF')
        plt.plot(cv_czar, pmf_total, color = 'tab:red', label = 'eABF + GaMDD')
        plt.plot(cv_czar, pmf_czar, color = 'tab:blue', label = 'eABF part')
        plt.plot(cv_czar, pmf_amd, color = 'tab:green', label = 'GaMDD part')
        plt.xlabel('Distance (nm)')
        plt.ylabel('$\Delta G$ (kcal/mol)')
        ax = plt.gca()
        ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
        ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
        ax.xaxis.get_major_formatter()._usetex = False
        ax.yaxis.get_major_formatter()._usetex = False
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim(0, 30)
        plt.legend(prop = {'size': 14}, fancybox = False, frameon = False)
        plt.tight_layout(pad = 0.2)
        output_png = outputname + '.png'
        plt.savefig(output_png, dpi=600, transparent=False)
        plt.close()

def get_total_count(countname):
    data = np.genfromtxt(countname, unpack = True)
    total_count = np.sum(data[-1])
    print(f'Total count for {countname}: {total_count:.2f}')
    return total_count

mergepmf('free_energy+1.production_abf.czar.pmf', 'free_energy+1.production_amd.cumulant.pmf', 'merge2', False)
#mergepmf('free_energy.production_abf.czar.pmf', 'free_energy.production_amd.reweight.pmf', 'merge1', False)
#get_total_count('free_energy.production_abf.zcount')
#get_total_count('free_energy.production_amd.count')
#get_total_count('deca.abf1.zcount')

#for i in range(2, 302):
    #eabf_pmf_file = 'deca.production.hist.czar_' + str(i).zfill(4) + '.pmf'
    #gamd_cumulant_pmf = 'deca.reweightamd1.cumulant.pmf_' + str(i + 300).zfill(4) + '.hist'
    #merge_output = 'merge_' + str(i).zfill(4)
    #mergepmf(eabf_pmf_file, gamd_cumulant_pmf, merge_output, plot = False)

#c1 = get_total_count('deca.production.hist_0000.zcount')
#c1 = get_total_count('deca.production.hist_0255.zcount')
#c2 = get_total_count('deca.production.hist_0256.zcount')
#c2 = get_total_count('deca.production.hist_0599.zcount')
#c3 = get_total_count('deca.abf1.zcount')
#print(c2 + c1)

#for i in range(0, 600):
    #count_name = 'deca.production.hist_' + str(i).zfill(4) + '.zcount'
    #c = get_total_count(count_name);
    #print(f'{i} {c:.2f}')

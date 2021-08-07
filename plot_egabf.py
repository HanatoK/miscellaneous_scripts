#!/usr/bin/env python3
from numpy import ceil
import pandas as pd
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
from matplotlib.ticker import AutoMinorLocator, AutoLocator

def plot_1d(data, outputname, xlabel, ylabel):
    plt.rcParams.update({
        "pgf.texsystem": "lualatex",
        "font.family": "serif",  # use serif/main font for text elements
        "text.usetex": True,     # use inline math for ticks
        "pgf.rcfonts": False,    # don't setup fonts from rc parameters
        "axes.labelsize": 28,
        "axes.linewidth": 2.0,
        'axes.unicode_minus': False,
        "font.size": 24,
        "pgf.preamble": '\n'.join([
            "\\usepackage{units}",
            "\\usepackage{metalogo}",
            "\\usepackage{unicode-math}",
            r"\setmathfont{MathJax_Math}",
            r"\setmainfont{FreeSans}",
        ])
    })
    w, h = figaspect(1/1.1)
    plt.figure(figsize = (w,h))
    plt.plot(data[0], data[1])
    # plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    ax = plt.gca()
    ax.set_xlim(-180.0, 180.0)
    ax.set_ylim(0, 16.0)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.tick_params(direction='in', which='major', length=6.0,
                   width=1.0, top=True, right=True, pad=8.0)
    ax.tick_params(direction='in', which='minor', length=3.0,
                   width=1.0, top=True, right=True, pad=8.0)
    ax.xaxis.set_major_locator(plt.MultipleLocator(90))
    ax.yaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.savefig(outputname, dpi=300, bbox_inches='tight', transparent=False)
    plt.close()


if __name__ == '__main__':
    phi1_data = pd.read_csv('trialanine_fes.abf1.czar.pmf',
                            delimiter='\s+', header=None, comment='#')
    phi2_data = pd.read_csv('trialanine_fes.abf3.czar.pmf',
                            delimiter='\s+', header=None, comment='#')
    phi3_data = pd.read_csv('trialanine_fes.abf5.czar.pmf',
                            delimiter='\s+', header=None, comment='#')
    psi1_data = pd.read_csv('trialanine_fes.abf2.czar.pmf',
                            delimiter='\s+', header=None, comment='#')
    psi2_data = pd.read_csv('trialanine_fes.abf4.czar.pmf',
                            delimiter='\s+', header=None, comment='#')
    psi3_data = pd.read_csv('trialanine_fes.abf6.czar.pmf',
                            delimiter='\s+', header=None, comment='#')
    plot_1d(data=phi1_data, outputname='PMF_phi1.png', xlabel=r'$\phi_1$ (°)', ylabel=r'$\Delta G$ (kcal/mol)')
    plot_1d(data=phi2_data, outputname='PMF_phi2.png', xlabel=r'$\phi_2$ (°)', ylabel=r'$\Delta G$ (kcal/mol)')
    plot_1d(data=phi3_data, outputname='PMF_phi3.png', xlabel=r'$\phi_3$ (°)', ylabel=r'$\Delta G$ (kcal/mol)')
    plot_1d(data=psi1_data, outputname='PMF_psi1.png', xlabel=r'$\psi_1$ (°)', ylabel=r'$\Delta G$ (kcal/mol)')
    plot_1d(data=psi2_data, outputname='PMF_psi2.png', xlabel=r'$\psi_2$ (°)', ylabel=r'$\Delta G$ (kcal/mol)')
    plot_1d(data=psi3_data, outputname='PMF_psi3.png', xlabel=r'$\psi_3$ (°)', ylabel=r'$\Delta G$ (kcal/mol)')

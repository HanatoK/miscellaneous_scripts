#!/usr/bin/env python3
import matplotlib
import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.use("pgf")
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 32,
    "axes.linewidth": 2.0,
    "font.size": 28,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{Arimo}",  # serif font via preamble
         ])
})
    

def read_2d_pmf(pmf_filename):
    pmf = pd.read_csv(
        pmf_filename, delimiter=r'\s+', header=None, comment='#')
    x = pmf[0].to_numpy()
    y = pmf[1].to_numpy()
    z = pmf[2].to_numpy()
    z = np.clip(z, 0, 16.0)
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    return xi, yi, zi


if __name__ == '__main__':
    minima_traj = [pd.read_csv(f'minima_{i}.csv') for i in range(1, 1+5)]
    xi, yi, zi = read_2d_pmf('alad.abf1.czar.pmf')
    xr, yr, zr = read_2d_pmf('alad_300ns.czar.pmf')
    fig = plt.figure(figsize=(24, 12))
    ax1 = fig.add_subplot(121)
    ax1.xaxis.get_major_formatter()._usetex = False
    ax1.yaxis.get_major_formatter()._usetex = False
    ax1.set_xlabel('CV1')
    ax1.set_ylabel('CV2')
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(6))
    color_map = plt.cm.get_cmap('tab20')
    cf1 = ax1.contourf(xi, yi, zi, 25, cmap='viridis')
    scatters1 = [ax1.scatter(m['CV1'].values, m['CV2'].values,
                             s=1.0, color=color_map(i))
                 for i, m in enumerate(minima_traj, start=1)]
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb1 = fig.colorbar(cf1, cax=cax, orientation='vertical')
    cb1.ax.xaxis.get_major_formatter()._usetex = False
    cb1.ax.yaxis.get_major_formatter()._usetex = False
    # second subplot
    ax2 = fig.add_subplot(122)
    cf2 = ax2.contourf(xr, yr, zr, 25, cmap='viridis')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb2 = fig.colorbar(cf2, cax=cax, orientation='vertical')
    cb2.ax.xaxis.get_major_formatter()._usetex = False
    cb2.ax.yaxis.get_major_formatter()._usetex = False
    scatters2 = [ax2.scatter(np.degrees(m['rad_phi'].values),
                             np.degrees(m['rad_psi'].values),
                             s=1.0, color=color_map(i))
                for i, m in enumerate(minima_traj, start=1)]
    ax2.xaxis.get_major_formatter()._usetex = False
    ax2.yaxis.get_major_formatter()._usetex = False
    ax2.set_xlim(-180.0, 180.0)
    ax2.set_ylim(-180.0, 180.0)
    ax2.set_xlabel(r'$\phi$')
    ax2.set_ylabel(r'$\psi$')
    ax2.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(6))
    plt.subplots_adjust(wspace=0.3)
    plt.savefig('test.png', dpi=300, bbox_inches='tight')

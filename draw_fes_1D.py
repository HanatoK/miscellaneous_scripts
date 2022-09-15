#!/usr/bin/env python3
import matplotlib
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import figaspect
from matplotlib.ticker import AutoMinorLocator
matplotlib.use("Agg")


def load_matplotlib_config():
    plt.rcParams.update({
        "pgf.texsystem": "lualatex",
        "font.family": "FreeSans",  # use serif/main font for text elements
        "mathtext.fontset": "stix",
        "text.usetex": False,     # use inline math for ticks
        "pgf.rcfonts": False,    # don't setup fonts from rc parameters
        "axes.labelsize": 28,
        "axes.linewidth": 2.0,
        "font.size": 24,
        'axes.unicode_minus': False
    })


def plot_1d_pmf(pmf_filename, outputname):
    data = pd.read_csv(pmf_filename, delimiter=r"\s+", comment="#", header=None)
    x = data[0].to_numpy()
    y = data[1].to_numpy()
    y = y - np.min(y)
    w, h = figaspect(1/1.1)
    plt.figure(figsize=(w, h))
    plt.plot(x, y)
    plt.xlabel(r"$\xi$")
    plt.ylabel(r"$\Delta G$ (kcal/mol)")
    ax = plt.gca()
    ax.tick_params(direction='in', which='major', length=4.0, width=1.0, top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=2.0, width=1.0, top=True, right=True)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.savefig(outputname, dpi=400, bbox_inches='tight', transparent=False)
    return y


def main():
    parser = argparse.ArgumentParser("Plot 1D PMF file")
    parser.add_argument("pmf", help="specify the PMF file")
    parser.add_argument("-o", "--output", help="specify the PNG output image file")
    args = parser.parse_args()
    if args.output is None:
        outputname = args.pmf + ".png"
    else:
        outputname = args.output
    load_matplotlib_config()
    plot_1d_pmf(args.pmf, outputname)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
import matplotlib
import os
import argparse
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "FreeSans",  # use serif/main font for text elements
    "text.usetex": False,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "mathtext.fontset": "stix",
    "axes.labelsize": 24,
    "axes.linewidth": 2.0,
    "font.size": 20,
    "axes.unicode_minus": False,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",
    ])
})
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker

class Axis:

    def __init__(self, lower_bound=0.0, upper_bound=0.0, bins=0, periodic=False):
        self.lowerBound = lower_bound
        self.upperBound = upper_bound
        self.bins = bins
        self.periodic = periodic
        if bins != 0:
            self.width = (self.upperBound - self.lowerBound) / float(self.bins)
        else:
            self.width = 0
        self.periodicLowerBound = self.lowerBound
        self.periodicUpperBound = self.upperBound

    def set_periodicity(self, periodic, periodic_lower, periodic_upper):
        self.periodic = periodic
        self.periodicUpperBound = periodic_upper
        self.periodicLowerBound = periodic_lower

    def get_lower_bound(self):
        return self.lowerBound

    def set_lower_bound(self, new_lower_bound):
        if self.width > 0:
            self.bins = round((self.upperBound - new_lower_bound) / self.width)
        else:
            self.bins = 0
        if self.bins == 0:
            self.lowerBound = new_lower_bound
        else:
            self.lowerBound = self.upperBound - float(self.bins) * self.width

    def get_upper_bound(self):
        return self.upperBound

    def set_upper_bound(self, new_upper_bound):
        if self.width > 0:
            self.bins = round((new_upper_bound - self.lowerBound) / self.width)
        else:
            self.bins = 0
        if self.bins == 0:
            self.upperBound = new_upper_bound
        else:
            self.upperBound = self.lowerBound + float(self.bins) * self.width

    def get_width(self):
        return self.width

    def set_width(self, new_width):
        self.bins = int(round((self.upperBound - self.lowerBound) / new_width))
        if self.bins == 0:
            self.width = new_width
        else:
            self.width = (self.upperBound - self.lowerBound) / float(self.bins)

    def get_bin(self):
        return self.bins

    def set_bin(self, new_bin):
        self.bins = new_bin
        self.width = (self.upperBound - self.lowerBound) / float(self.bins)

    def in_boundary(self, x):
        x = self.wrap(x)
        if (x < self.lowerBound) or (x > self.upperBound):
            return False
        else:
            return True

    def wrap(self, x):
        from math import fabs, isclose
        if not self.periodic:
            return x
        if (x >= self.periodicLowerBound) and (x <= self.periodicUpperBound):
            return x
        periodicity = self.periodicUpperBound - self.periodicLowerBound
        if x < self.periodicLowerBound:
            dist_to_lower = self.periodicLowerBound - x
            num_period_add = int(dist_to_lower / periodicity)
            tmp = fabs(dist_to_lower / periodicity - round(dist_to_lower / periodicity))
            if isclose(tmp, 0.0):
                x += num_period_add * periodicity
            else:
                x += (num_period_add + 1) * periodicity
        if x > self.periodicUpperBound:
            dist_to_upper = x - self.periodicUpperBound
            num_period_subtract = int(dist_to_upper / periodicity)
            tmp = fabs(dist_to_upper / periodicity - round(dist_to_upper / periodicity))
            if isclose(tmp, 0.0):
                x -= num_period_subtract * periodicity
            else:
                x -= (num_period_subtract + 1) * periodicity
        return x

    def index(self, x, boundary_check=False):
        from math import floor
        x = self.wrap(x)
        check_result = True
        if boundary_check is True:
            check_result = self.in_boundary(x)
        if check_result is False:
            return 0, check_result
        idx = int(floor((x - self.lowerBound) / self.width))
        if idx == self.bins:
            idx = idx - 1
        return idx, check_result

    def dist(self, x, reference):
        from math import fabs
        if not self.periodic:
            return x - reference
        else:
            x = self.wrap(x)
            reference = self.wrap(reference)
            dist = x - reference
            p = self.period()
            if fabs(dist) > (p * 0.5):
                if reference > x:
                    return dist + p
                elif reference < x:
                    return dist - p
                else:
                    return dist
            else:
                return dist

    def period(self):
        return self.periodicUpperBound - self.periodicLowerBound

    def get_middle_points(self):
        result = [self.lowerBound + (i + 0.5) * self.width for i in range(0, self.bins)]
        return result

    def info_header(self):
        pbc = 0
        if self.periodic is True:
            pbc = 1
        return f'# {self.lowerBound:.9f} {self.width:.9f} {self.bins:d} {pbc}'

    def __str__(self) -> str:
        s = f'boundary: [{self.lowerBound}, {self.upperBound}] ; width: {self.width} ; PBC: {self.periodic}'
        return s
    
    @staticmethod
    def from_str(s):
        tmp = s.split()
        if len(tmp) < 5:
            raise RuntimeError('Incorrect format of axis.')
        lower_bound = float(tmp[1])
        width = float(tmp[2])
        num_bins = int(tmp[3])
        upper_bound = float(lower_bound + width * num_bins)
        pbc = True
        if int(tmp[4]) == 0:
            pbc = False
        ax = Axis(lower_bound=lower_bound, upper_bound=upper_bound,
                  bins=num_bins, periodic=pbc)
        if pbc is True:
            ax.set_periodicity(periodic=pbc, periodic_lower=lower_bound,
                               periodic_upper=upper_bound)
        return ax


def parse_header(gradfile):
    # read header
    with open(gradfile, 'r') as f_input:
        dim = int(f_input.readline().split()[1])
        if dim != 1:
            raise RuntimeError('Not a valid 1D grad file.')
        ax = Axis.from_str(f_input.readline())
        # since this is 1D integration by the trapezoidal rule, increasing 1 bin
        ax.lowerBound = ax.lowerBound - 0.5 * ax.get_width()
        ax.upperBound = ax.upperBound + 0.5 * ax.get_width()
        ax.bins += 1
        return ax


def integrate1D(x, grad):
    width = x[1] - x[0]
    start_x = x[0] - 0.5 * width
    start_pmf = 0
    pmf_cv = [start_x]
    pmf = [start_pmf]
    for g_i in np.nditer(grad):
        start_x += width
        start_pmf += g_i * width
        pmf_cv.append(start_x)
        pmf.append(start_pmf)
    return np.array(pmf_cv), np.array(pmf)


def plot_fes_1D(cv, pmf, xtitle, ytitle, label, color, outputname):
    w, h = figaspect(1/1.1)
    plt.figure(figsize = (w,h))
    alpha_ratio = 1.0
    plt.plot(cv, pmf, label = label, color = color, alpha = alpha_ratio)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    ax = plt.gca()
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
    ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
    plt.legend(prop = {'size': 14}, fancybox = False, frameon = False)
    plt.tight_layout()
    plt.savefig(outputname, dpi=600, transparent=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("grad", help = "specify the gradient file")
    parser.add_argument("-o", "--output", help = "specify the output prefix")
    args = parser.parse_args()
    header = parse_header(args.grad)
    x, grad_x = np.genfromtxt(args.grad, unpack = True)
    cv, pmf = integrate1D(x, grad_x)
    pmf = pmf - np.min(pmf)
    outputprefix = args.output
    outputpmf = outputprefix + '.pmf'
    outputpng = outputprefix + '.png'
    with open(outputpmf, 'w') as f_pmf:
        f_pmf.write('# 1\n')
        f_pmf.write(header.info_header() + '\n')
        np.savetxt(f_pmf, np.c_[cv, pmf], fmt = '%10.2f %12.7f')
    plot_fes_1D(cv, pmf, r'$\xi$', '$\Delta G$ (kcal/mol)', None, 'red', outputpng)


if __name__ == '__main__':
    main()

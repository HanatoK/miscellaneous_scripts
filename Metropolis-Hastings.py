#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 20,
    "axes.linewidth": 2.0,
    'axes.unicode_minus': False,
    "font.size": 16,
    "pgf.preamble": '\n'.join([
        "\\usepackage{units}",
        "\\usepackage{metalogo}",
        "\\usepackage{unicode-math}",
        r"\setmathfont{MathJax_Math}",
        r"\setmainfont{FreeSans}",
    ])
})

def target_function(x):
    if x < 0:
        return 0
    else:
        return np.exp(-1.0 * x)

def adjacent_average(X):
    Y = list()
    if len(X) <= 1:
        return X
    else:
        for i in range(1, len(X)):
            Y.append(0.5 * (X[i] + X[i-1]))
    return np.array(Y)

if __name__ == '__main__':
    # save to a list
    traj = list()
    # initial value
    x = np.random.rand()
    # number of steps
    N = 200000
    # running
    current_x = x
    traj.append(current_x)
    for i in range(N):
        # random drift
        drift = np.random.normal()
        # proposed new
        proposed_x = current_x + drift
        # calculate the probability
        p_proposed = target_function(proposed_x)
        p_current = target_function(current_x)
        # calculate the acceptance ratio
        acceptance_ratio = min(1.0, p_proposed / p_current)
        if np.random.uniform(low=0.0, high=1.0) < acceptance_ratio:
            # accept
            current_x = proposed_x
        else:
            # reject
            pass
        traj.append(current_x)
    # save to file
    with open('traj.dat', 'w') as fOutput:
        for x in traj:
            fOutput.write(f'{x:12.7f}\n')
    # plotting
    data, edges = np.histogram(traj, bins=100, range=(0, 10), density=True)
    bin_centers = adjacent_average(edges)
    plt.plot(bin_centers, data, label='Sampling')
    plt.plot(bin_centers, np.exp(-1.0 * bin_centers), label='Actual\ndistribution')
    ax = plt.gca()
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.tick_params(direction='in', which='major', length=6.0,
                    width=1.0, top=True, right=True, pad=8.0)
    ax.tick_params(direction='in', which='minor', length=3.0,
                    width=1.0, top=True, right=True, pad=8.0)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(0, 10)
    plt.legend(prop = {'size': 16}, fancybox = False, frameon = False, handlelength=1.0, handletextpad=0.5)
    plt.savefig('test.png', dpi=300, bbox_inches='tight')
    np.savetxt('hist.dat', np.c_[data, bin_centers], fmt='%12.7f')
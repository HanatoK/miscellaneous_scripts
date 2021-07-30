#!/usr/bin/env python3
# import numpy as np
# import pandas as pd
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",
    ])
})
from matplotlib.figure import figaspect

if __name__ == '__main__':
    data = []
    with open('committor_result.log', 'r') as fInput:
        for line in fInput:
            if line.startswith('gpath'):
                fields = line.split()
                data.append(float(fields[3]))
    w, h = figaspect(1/1.1)
    plt.figure(figsize = (w,h))
    x, bins, p = plt.hist(data, bins=20, range=(0.0, 1.0), edgecolor='black', linewidth=1.0, density=True)
    for item in p:
        item.set_height(item.get_height()/sum(x))
    plt.title('A-M1-M4-B')
    plt.xlabel('$p$($s(\mathbf{z})$ < 0.7)')
    plt.ylabel('Probability')
    ax = plt.gca()
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1)
    ax.tick_params(direction='in', which='major', length=6.0, width=1.0, top=True, right=True, pad=10.0)
    ax.tick_params(direction='in', which='minor', length=3.0, width=1.0, top=True, right=True, pad=10.0)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    plt.savefig('committor_result.png', dpi=300, transparent=False, bbox_inches='tight')
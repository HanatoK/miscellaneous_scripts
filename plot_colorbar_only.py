#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "FreeSans",  # use serif/main font for text elements
    "mathtext.fontset": "stix",
    "text.usetex": False,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
    'axes.unicode_minus': False,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",  # serif font via preamble
         ])
})


# see https://stackoverflow.com/questions/16595138/standalone-colorbar-matplotlib
img = plt.imshow(np.array([[0, 14]]), cmap='plasma_r')
plt.gca().set_visible(False)
# img.set_visible(False)

clb = plt.colorbar(orientation='horizontal')
clb.ax.xaxis.set_major_locator(plt.FixedLocator(np.linspace(0, 14, 8)))
plt.savefig('colorbar.png', dpi=300, bbox_inches='tight', transparent=True)

#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def hist_colvars_traj(inputfilename, outputfilename):
    colvars_traj = pd.read_csv(inputfilename, delimiter='\s+', header=None, comment='#')
    data = colvars_traj[4]
    plt.hist(data, bins=50, range=(0.489, 0.491), edgecolor='black', linewidth=0.5) 
    plt.savefig(outputfilename, dpi=300, transparent=False, bbox_inches='tight')
    return

if __name__ == '__main__':
    hist_colvars_traj('gpath_structures/frame_cv.dat', 'gpath_initial.png')

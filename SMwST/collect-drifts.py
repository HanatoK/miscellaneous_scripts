#!/usr/bin/env python3
from glob import glob
import argparse
import os
from shutil import copyfile
import numpy as np
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 36,
    "axes.linewidth": 2.0,
    "font.size": 28,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",  # serif font via preamble
         ])
})
from matplotlib.figure import figaspect
#from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('n', help='number of images')
parser.add_argument('--swarm_steps', type=int, default=10, help = 'steps of short unrestrained steps (lengths) of a swarm')
parser.add_argument('--equil_steps', type=int, default=20000, help = 'steps of the restrained equilibrium steps before forming a swarm')
args = parser.parse_args()
num_images = int(args.n)

# get the number of trajectories per swarm
output_dirs = glob('output/*')
num_copies = int(len(output_dirs) / num_images)
print(f'Number of images : {num_images}')
print(f'Number of trajectories per image : {num_copies}')
output_dirs.sort()
#result_dirs = output_dirs[::num_copies]

# SMwST parameters to determine the drift steps
num_drift_steps = args.swarm_steps
num_equilibrium_steps = args.equil_steps
timestep = 0.5 # fs

# helper function to extract drifts from a colvars trajectory
# return drifts of all iterations of a trajectory
def extract_drift(colvars_traj, num_drift_steps, num_equilibrium_steps):
    source = []
    target = []
    j = 0
    i = 1
    with open(colvars_traj, 'r') as finput:
        for line in finput:
            if line.startswith('#'):
                continue
            if line:
                fields = line.split()
                if ((int(fields[0]) - num_drift_steps * j) % num_equilibrium_steps == 0):
                    j = j + 1
                    target.append(np.array([float(fields[k]) for k in range(1, len(fields))]))
                if ((int(fields[0]) - num_drift_steps * i) % num_equilibrium_steps == 0):
                    i = i + 1
                    source.append(np.array([float(fields[k]) for k in range(1, len(fields))]))
    result = np.array(source) - np.array(target[:len(source)])
    result = np.sqrt(np.sum(result * result, axis=1))
    return result

# helper function for plotting
def plot_dSm(data, image_index, outputname):
    w, h = figaspect(1/2)
    plt.figure(figsize = (w,h))
    x = np.arange(0, len(data), 1)
    plt.plot(x, data)
    plt.title(f'm = {image_index+1}')
    plt.xlabel(r'Iterations')
    plt.ylabel(r'$\sqrt{S^m(\delta\tau)}/\delta\tau$ (°)')
    ax = plt.gca()
    ax.set_ylim(0, 1.5)
    ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
    ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    plt.savefig(outputname, dpi=300, bbox_inches = 'tight', transparent=False)
    plt.close()

# make directories and copy drifts
base_directory = os.getcwd() + '/drift_results'
if not os.path.exists(base_directory):
    os.makedirs(base_directory)

# search drifts from covlars traj and copy them
digits_images = len(str(num_images - 1))
digits_copies = len(str(num_copies - 1))
digits_simulations = len(str(num_images * num_copies - 1))
for i in range(0, num_images):
    image_dirname = base_directory + '/' + str(i).zfill(digits_images)
    if not os.path.exists(image_dirname):
        os.makedirs(image_dirname)
    # read trajectories in a swarm
    firsttime = True
    for j in range(0, num_copies):
        traj_filename = 'output/' + str(i * num_copies + j).zfill(digits_simulations) + '/trialanine.job0000.colvars.traj'
        drift_data = extract_drift(traj_filename, num_drift_steps, num_equilibrium_steps)
        outputname = image_dirname + '/' + str(j).zfill(digits_copies)
        np.savetxt(outputname, drift_data, fmt='%12.7f')
        # compute average
        if firsttime:
            sum_drift = drift_data
            sum_drift_square = drift_data * drift_data
            firsttime = False
        else:
            sum_drift += drift_data
    sum_drift /= num_copies
    sum_drift_square /= num_copies
    average_outputname = image_dirname + '/average.dat'
    square_average_outputname = image_dirname + '/square_average.dat'
    # this is what Benoit wants
    dSm_tau_t = np.sqrt(sum_drift_square) / (timestep * num_drift_steps)
    dSm_tau_t_outputname = image_dirname + '/dSm_tau_t.dat'
    dSm_tau_t_pngname = image_dirname + '/dSm_tau_t.png'
    np.savetxt(average_outputname, sum_drift, fmt='%12.7f')
    np.savetxt(square_average_outputname, sum_drift_square, fmt='%12.7f')
    np.savetxt(dSm_tau_t_outputname, dSm_tau_t, fmt='%12.7f')
    plot_dSm(dSm_tau_t, i, dSm_tau_t_pngname)

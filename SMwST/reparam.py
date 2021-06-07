#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "FreeSans",  # use serif/main font for text elements
    "text.usetex": False,    # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
})
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator

# calculate L(i)
def calc_l_i(image_coordinates, i):
    # image_coordinates is (M, N) numpy 2D array
    # where M is the number of images and N is the number of CVs
    # calculate |z^m - z^{m-1}|, the first element is 0
    if i == 0:
        return 0
    # calculate z^m - z^{m-1}
    dist_z = np.diff(image_coordinates[0:i+1], axis=0)
    # calculate norm
    dist_z = np.linalg.norm(dist_z, axis=1)
    # now accumulate sum_norm to l(i/M)
    L = np.sum(dist_z)
    return L


# a single raparametrization move
def reparam_move(image_coordinates):
    # first compute the total length
    total_L = calc_l_i(image_coordinates, len(image_coordinates)-1)
    M = float(len(image_coordinates)-1)
    # keep the first image fixed
    new_image_coordinates = [np.copy(image_coordinates[0])]
    # iterate over the middle images
    for i in range(1, len(image_coordinates)-1):
        vec = image_coordinates[i] - new_image_coordinates[i-1]
        # normalize it to a unit vector
        vec = vec / np.linalg.norm(vec)
        # new moving vector
        new_vec = (i / M * total_L - calc_l_i(new_image_coordinates, i-1)) * vec
        # append new image
        new_image_coordinates.append(new_image_coordinates[i-1] + new_vec)
    # append the fixed last image
    new_image_coordinates.append(np.copy(image_coordinates[-1]))
    return np.array(new_image_coordinates)


def path_rmsd(X, ref):
    diff2 = (X - ref) * (X - ref)
    rmsd = np.sqrt(np.sum(diff2) / len(X))
    return rmsd


def check_distances(image_coordinates):
    dist = np.linalg.norm(np.diff(image_coordinates, axis=0), axis=1)
    #print(dist)
    return np.std(dist)


# iteratively reparmetrization
def reparametrization(image_coordinates):
    total_L = calc_l_i(image_coordinates, len(image_coordinates)-1)
    max_iteration = 100
    previous_path = image_coordinates
    check_distances(previous_path)
    next_path = reparam_move(image_coordinates)
    iteration = 1
    while (iteration < 25):
        rmsd = path_rmsd(next_path, previous_path)
        stddev = check_distances(next_path)
        print(f'Iteration {iteration:03d}: RMSE = {rmsd:10.7f} ; length stddev = {stddev:10.7f} ; total length = {total_L:10.7f}')
        previous_path = next_path
        next_path = reparam_move(next_path)
        iteration = iteration + 1
        total_L = calc_l_i(next_path, len(image_coordinates)-1)
    return next_path


if __name__ == '__main__':
    N = 20
    X = np.linspace(-5, 5, N)
    ref = np.sin(2*X)
    Y = ref + np.append(np.append(0, np.random.normal(0, 0.3, N-2)), 0)
    coords = np.c_[X, Y]
    reparam_coords = reparametrization(coords)
    np.savetxt('origin.txt', coords, fmt='%12.7f')
    np.savetxt('reparam.txt', reparam_coords, fmt='%12.7f')
    # plotting
    w, h = figaspect(1/2)
    plt.plot(X, ref, label='sin(x)', color='red', alpha=0.5)
    plt.scatter(X, ref, marker='x', color='red')
    plt.plot(X, Y, label='Origin', color='orange', alpha=0.5)
    plt.scatter(X, Y, marker='x', color='orange')
    reparam_X = np.transpose(reparam_coords)[0]
    reparam_Y = np.transpose(reparam_coords)[1]
    plt.plot(reparam_X, reparam_Y, label='Reparametrization', color='green', alpha=0.5)
    plt.scatter(reparam_X, reparam_Y, marker='x', color='green')
    plt.xlabel('X')
    plt.ylabel('Y')
    ax = plt.gca()
    ax.tick_params(direction='in', which='major', length=6.0, width=1.0, top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=3.0, width=1.0, top=True, right=True)
    ax.set_ylim(-2.5, 2.5)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.legend(prop = {'size': 14}, fancybox = False, frameon = False)
    plt.savefig('reparam.png', dpi=300, bbox_inches='tight', transparent=False)

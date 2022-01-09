#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import CubicSpline
import argparse
import pandas as pd

def reparametrize(positions, num_images=None, resolution=1000):
    # do interpolation
    X = np.linspace(0, resolution, len(positions))
    print(X)
    fX = CubicSpline(X, positions, axis=0, bc_type='natural')
    if num_images is None:
        num_images = len(positions)
    new_X = np.linspace(0, resolution, num_images * resolution)
    new_Y = fX(new_X)
    # compute the total length
    adjacent_diff = np.diff(new_Y, axis=0)
    all_lengths = np.sqrt(np.sum(adjacent_diff * adjacent_diff, axis=1))
    total_length = np.sum(all_lengths)
    # find points on the spline
    adjacent_length = total_length / (num_images - 1)
    idx = 1
    sum_l = 0
    result = [new_Y[0]]
    for i, l in enumerate(all_lengths[:len(all_lengths)-1], start=1):
        sum_l += l
        if (sum_l >= idx * adjacent_length):
            result.append(new_Y[i])
            idx = idx + 1
    result.append(new_Y[-1])
    return np.array(result)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='path position file')
    parser.add_argument('-n', '--num_images', type=int, help='number of desired images (default to the same images of the path)')
    parser.add_argument('-i', '--num_iterations', default=1, type=int, help='number of iterations')
    parser.add_argument('-o', '--output', help='output filename')
    args = parser.parse_args()
    positions = pd.read_csv(args.path, delimiter='\s+', comment='#', header=None).to_numpy()
    Y = positions.copy()
    for i in range(0, args.num_iterations):
        Y = reparametrize(positions=Y, num_images=args.num_images)
    # write the output
    with open(args.output, 'w') as fOutput:
        for point in Y:
            for elem in point:
                fOutput.write(f' {elem:15.10f}')
            fOutput.write('\n')
    # debug: show neighboring distances
    adjacent_diff = np.diff(Y, axis=0)
    all_lengths = np.sqrt(np.sum(adjacent_diff * adjacent_diff, axis=1))
    for i, l in enumerate(all_lengths):
        print(f"Distance between image {i} and {i+1}: {l:12.5f}")
    # debug: show moved distances
    dist_Y_pos = Y - positions
    distances = np.sqrt(np.sum(dist_Y_pos * dist_Y_pos, axis=1))
    for i, l in enumerate(distances):
        print(f'Image {i} has moved: {l:12.5f}')
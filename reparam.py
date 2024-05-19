#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import CubicSpline
import argparse
import pandas as pd


def reparametrize(positions, num_images=None, resolution=1000):
    # do interpolation
    X = np.linspace(0, resolution, len(positions))
    # print(X)
    fX = CubicSpline(X, positions, axis=0, bc_type='natural')
    df = fX.derivative()
    if num_images is None:
        num_images = len(positions)
    new_X = np.linspace(0, resolution, num_images * resolution)
    new_Y = fX(new_X)
    new_dY = df(new_X)
    # compute the total length
    adjacent_diff = np.diff(new_Y, axis=0)
    all_lengths = np.sqrt(np.sum(adjacent_diff * adjacent_diff, axis=1))
    total_length = np.sum(all_lengths)
    # find points on the spline
    adjacent_length = total_length / (num_images - 1)
    cumsum_l = np.cumsum(all_lengths[:-1])
    cumsum_adj_len = np.arange(1, num_images - 1) * adjacent_length
    idx = np.searchsorted(cumsum_l, cumsum_adj_len) + 1
    result = np.concatenate((np.asarray([new_Y[0]]), new_Y[idx], np.asarray([new_Y[-1]])))
    derivative = np.concatenate((np.asarray([new_dY[0]]), new_dY[idx], np.asarray([new_dY[-1]])))
    return np.asarray(result), np.asarray(derivative)


# translate from C++
def remove_loops(positions):
    def image_distance(pos_x, pos_y):
        diff = np.asarray(pos_x) - np.asarray(pos_y)
        return np.sqrt(np.sum(diff * diff))
    results = list(positions)
    has_loop = False
    while True:
        has_loop = False
        for i in range(1, len(results) - 1):
            neighbor_dist = image_distance(results[i-1], results[i])
            j = i + 1
            while j < len(results):
                dist = image_distance(results[i-1], results[j])
                if dist < neighbor_dist:
                    has_loop = True
                    break
                j = j + 1
            if has_loop:
                print(f'Remove indexes from {i} to {j}')
                del results[i:j]
                break
        if not has_loop:
            break
    return np.array(results)


def check_equidistance(positions):
    adjacent_diff = np.diff(positions, axis=0)
    all_lengths = np.sqrt(np.sum(adjacent_diff * adjacent_diff, axis=1))
    for i, l in enumerate(all_lengths):
        print(f"Distance between image {i} and {i+1}: {l:12.5f}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='path position file')
    parser.add_argument('-n', '--num_images', type=int, help='number of desired images (default to the same images of the path)')
    parser.add_argument('-i', '--num_iterations', default=1, type=int, help='number of iterations')
    parser.add_argument('-o', '--output', help='output filename')
    args = parser.parse_args()
    positions = pd.read_csv(args.path, delimiter=r'\s+', comment='#', header=None).to_numpy()
    Y = positions.copy()
    for i in range(0, args.num_iterations):
        Y = remove_loops(Y)
        Y, dY = reparametrize(positions=Y, num_images=args.num_images)
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
    if np.shape(Y) == np.shape(positions):
        dist_Y_pos = Y - positions
        distances = np.sqrt(np.sum(dist_Y_pos * dist_Y_pos, axis=1))
        for i, l in enumerate(distances):
            print(f'Image {i} has moved: {l:12.5f}')
    # debug: show the derivatives
    print('Derivative from cubic splines:')
    # np.set_printoptions(formatter={'float': '{: 12.5ff}'.format})
    print(dY)


if __name__ == '__main__':
    main()

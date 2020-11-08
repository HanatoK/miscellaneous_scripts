#!/usr/bin/env python3
import math
import numpy as np


def diffSquare(x1, x2, isPeriodic, period):
    if isPeriodic is False:
        return (x1 - x2) * (x1 - x2)
    else:
        diff = max(x1, x2) - min(x1, x2)
        halfPeriod = 0.5 * period
        if diff > halfPeriod:
            return (diff - halfPeriod) * (diff - halfPeriod)
        else:
            return diff * diff


def nodeDistance(l1, l2):
    dist = 0
    for (a, b) in zip(l1, l2):
        dist = dist + diffSquare(a, b, True, 360.0)
    dist = math.sqrt(dist)
    return dist


num_iterations = 100
reference_iteration = f'path_{num_iterations-1:03d}'
for i in range(1, num_iterations):
    previous_iteration = np.genfromtxt(f'path_{i-1:03d}', unpack=True).transpose()
    current_iteration = np.genfromtxt(reference_iteration, unpack=True).transpose()
    rmsd = 0
    for image_previous, image_current in zip(previous_iteration, current_iteration):
        dist = nodeDistance(image_previous, image_current)
        rmsd += dist * dist
    rmsd /= len(current_iteration)
    rmsd = math.sqrt(rmsd)
    print(f'RMSD between {i-1:3d} and {reference_iteration} is {rmsd:12.5f}')

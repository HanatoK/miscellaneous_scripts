#!/usr/bin/env python3
import numpy as np
import argparse

def images_rmsd(data):
    num_images = len(data)
    rmsd = []
    for i in range(1, num_images):
        diff = data[i] - data[i-1]
        dist = np.sqrt(np.dot(diff, diff))
        print(f'Distance between image[{i-1:03d}] and image[{i:03d}] = {dist:12.7f}')
        rmsd.append(dist)
    return np.array(rmsd)

parser = argparse.ArgumentParser()
parser.add_argument('path', nargs = '+', help = 'path file(s)')
args = parser.parse_args()
path = args.path

for p in path:
    data = np.genfromtxt(p, unpack=True)
    images_rmsd(data.transpose())

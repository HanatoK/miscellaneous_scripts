#!/usr/bin/env python3
import numpy as np
import os
import argparse

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('count', nargs = '+', help = 'specify the count file(s)')
parser.add_argument('-t', '--timestep', help = 'timestep (fs) for compute the estimated simulation time')
args = parser.parse_args()

def sum_count(count_file, timestep = 2.0):
    data = np.genfromtxt(count_file, unpack = True)
    total_count = np.sum(data[-1])
    estimated_time = timestep * total_count / 1e6
    print(f'Total count of {count_file} is {total_count:.2f} ; Estimated simulation time : {estimated_time:.2f} (ns)')
    return total_count

for infile in args.count:
    if args.timestep:
        sum_count(infile, timestep = float(args.timestep))
    else:
        sum_count(infile)

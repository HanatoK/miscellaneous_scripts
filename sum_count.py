#!/usr/bin/env python3
# import numpy as np
# import os
import argparse

def sum_count(count_file):
    total_count = 0
    max_count = 0
    first_time = True
    with open(count_file, 'r') as f_input:
        for line in f_input:
            if line.startswith('#'):
                continue
            fields = line.split()
            if fields:
                if first_time:
                    print(f'The count file {count_file} has {len(fields)} data columns.')
                    first_time = False
                c = float(fields[-1])
                total_count += c
                if c > max_count:
                    max_count = c
    return total_count, max_count

if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('count', nargs = '+', help = 'specify the count file(s)')
    parser.add_argument('-t', '--timestep', type=float, default=2.0,
                        help='timestep (fs) for compute the estimated simulation time')
    args = parser.parse_args()
    for infile in args.count:
        if args.timestep:
            total_count, max_count = sum_count(infile)
            estimated_time = args.timestep * total_count / 1e6
            print(f'Total count of {infile} is {total_count:.2f}\n'
                  f'Estimated simulation time : {estimated_time:.2f} (ns)\n'
                  f'Max count = {max_count}')
        else:
            total_count, max_count = sum_count(infile)
            print(f'Total count of {infile} is {total_count:.2f}\n'
                  f'Max count = {max_count}')

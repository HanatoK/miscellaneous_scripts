#!/usr/bin/env python3
import os
import argparse

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--prefix', default='output/', help='the prefix of directories')
parser.add_argument('n', help='number of directories')
args = parser.parse_args()
num_images = int(args.n)
dir_prefix = args.prefix

# get the number of digits
num_digits = len(str(num_images - 1))
# make directories
for i in range(0, num_images):
    output_dir = dir_prefix + str(i).zfill(num_digits)
    print(f'--> Creating directory: {output_dir}')
    os.makedirs(output_dir)

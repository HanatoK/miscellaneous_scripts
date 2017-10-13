#!/usr/bin/python3
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('data', nargs = '+', help = 'specify the data file(s)')
parser.add_argument('-o', '--output', help = 'PNG output file')
parser.add_argument('-n', '--column', help = 'column you want to calculate')
args = parser.parse_args()

datafile = args.data
col = int(args.column) - 1
outfile = args.output
#out = open(outfile, 'w')
first = True
for infile in datafile:
    if first is True:
        x = np.genfromtxt(infile, unpack=True, usecols=col)
        first = False
    x = x + np.genfromtxt(infile, unpack=True, usecols=col)
x = x / len(datafile)
idx = np.linspace(1, len(x), len(x)).astype(int)
np.savetxt(outfile, np.column_stack([idx, x]), fmt=['%d','%10.7f'], delimiter='\t') 

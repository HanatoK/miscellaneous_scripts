#!/usr/bin/python3
import numpy as np
import argparse

def rmsd(ref, x):
    deviation = np.array(x) - np.array(ref)
    rmsd = np.sqrt(np.mean(deviation * deviation))
    return rmsd

parser = argparse.ArgumentParser()
parser.add_argument('pmf', nargs = '+', help = 'specify the pmf file(s)')
parser.add_argument('-r', '--reference', help = 'specify the reference PMF')
parser.add_argument('-n', '--column', help = 'column you want to calculate')
parser.add_argument('-o', '--output', help = 'RMSD data output file')
args = parser.parse_args()
ref = args.reference
pmf = args.pmf
col = int(args.column) - 1
if args.output is None:
    out = 'fesrmsd.dat'
else:
    out = args.output
refcol = np.genfromtxt(ref, unpack=True, usecols=col)
outputfile = open(out, 'w')
for file in pmf:
    x = np.genfromtxt(file, unpack=True, usecols=col)
    line = file + ' ' + str(rmsd(refcol, x)) + '\n'
    outputfile.write(line)
    print(line, end='')
outputfile.close()
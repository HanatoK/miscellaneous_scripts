#!/usr/bin/env python3
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('pmf', nargs = '+', help = 'specify the pmf file(s)')
args = parser.parse_args()
pmffiles = args.pmf
firsttime = True
count = 0
for fn in pmffiles:
    data = np.genfromtxt(fn, unpack=True)
    lastcol = len(data) - 1
    if firsttime is True:
        sumpmf = data[lastcol]
        count = count + 1
        firsttime = False
    else:
        sumpmf += data[lastcol]
        count = count + 1
    avgpmf = sumpmf / count
    outfile = "new_" + os.path.basename(fn)
    data[lastcol] = avgpmf
    np.savetxt(outfile, np.transpose(data), fmt='%10.7f')

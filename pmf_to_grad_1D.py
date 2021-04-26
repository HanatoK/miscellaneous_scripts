#!/usr/bin/env python3
import numpy as np
import os
import argparse

def pmf_to_grad_1D(pmf_filename, grad_filename):
    first_dataline = True
    previous_x = None
    previous_pmf = None
    with open(pmf_filename, 'r') as fInput, open(grad_filename, 'w') as fOutput:
        for line in fInput:
            if line.startswith('#'):
                fOutput.write(line)
            else:
                fields = line.split()
                if fields:
                    if first_dataline:
                        previous_x = float(fields[0])
                        previous_pmf = float(fields[1])
                        first_dataline = False
                    else:
                        current_x = float(fields[0])
                        current_pmf = float(fields[1])
                        width = current_x - previous_x
                        grad = (current_pmf - previous_pmf) / width
                        grad_x = 0.5 * (current_x + previous_x)
                        fOutput.write(f'{grad_x:12.7f} {grad:12.7f}\n')
                        previous_x = current_x
                        previous_pmf = current_pmf
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pmf', help='PMF file')
    parser.add_argument('-o', '--output', default='output.grad', help='output gradient file')
    args = parser.parse_args()
    pmf_to_grad_1D(args.pmf, args.output)

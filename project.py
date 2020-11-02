#!/usr/bin/env python3
import numpy as np
import argparse


def project(PCV_vectors_filename, input_filename):
    source_data = np.transpose(np.genfromtxt(input_filename, unpack=True))
    PCV_vectors = np.transpose(np.genfromtxt(PCV_vectors_filename, unpack=True))
    return np.dot(source_data, PCV_vectors)


def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('vectors', help='PCV vectors file')
    parser.add_argument('-i', '--input', help='the input data')
    parser.add_argument('-o', '--output', default='output', help='the output file')
    args = parser.parse_args()
    PCV_vectors_filename = args.vectors
    input_filename = args.input
    output_filename = args.output

    projected_data = project(PCV_vectors_filename, input_filename)
    np.savetxt(output_filename, projected_data, fmt='%12.7f')


if __name__ == "__main__":
    main()

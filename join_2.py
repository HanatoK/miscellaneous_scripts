#!/usr/bin/env python3
import glob
import argparse
from natsort import natsorted, ns


def join_files(input_filename_list, output_filename):
    print(f'Will join {len(input_filename_list)} files')
    with open(output_filename, 'w') as output_handle:
        for input_filename in input_filename_list:
            print(f'Reading {input_filename}')
            with open(input_filename, 'r') as input_handle:
                for line in input_handle:
                    output_handle.write(line)


def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help = 'input files')
    parser.add_argument("-o", "--output", default = "output", help = "output file")
    args = parser.parse_args()
    input_files = args.files
    output_filename = args.output

    # sort the list
    input_filename_list = natsorted(input_files)
    join_files(input_filename_list, output_filename)


if __name__ == "__main__":
    main()

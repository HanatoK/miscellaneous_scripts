#!/usr/bin/env python3
import glob
import argparse


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
    parser.add_argument('prefix', help = 'input prefix')
    parser.add_argument("-o", "--output", default = "output", help = "the output file")
    args = parser.parse_args()
    prefix = args.prefix
    output_filename = args.output

    # globbing
    input_filename_list = sorted(glob.glob(prefix+'*'))
    join_files(input_filename_list, output_filename)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse


def detect_rows(history_file):
    first_time = True
    line_count = 0
    with open(history_file, 'r') as f_input:
        for line in f_input:
            line = line.strip()
            if line.startswith('#'):
                if first_time is False:
                    break
                continue
            else:
                if line:
                    first_time = False
                    line_count += 1
    return line_count


def read_history_file(history_file, num_rows_to_split):
    data = pd.read_csv(history_file, delimiter=r'\s+', comment='#', header=None)
    # check if the number of total history data lines are multiples of target lines
    if data.shape[0] % num_rows_to_split != 0:
        raise RuntimeError(f'Number of history lines ({data.shape[0]}) is '
                           f'not multiple of target lines ({num_rows_to_split})!')
    # split the data
    list_of_data = [data.iloc[start:start+num_rows_to_split] for start in range(0, data.shape[0], num_rows_to_split)]
    return list_of_data


def compute_rmsd(source_grad, reference_grad):
    # the gradients should be in the last d/2 cols of the table
    grad_col_start = source_grad.shape[1] // 2
    X = source_grad.iloc[:, grad_col_start:].to_numpy()
    X_ref = reference_grad.iloc[:, grad_col_start:].to_numpy()
    diff2 = (X - X_ref) ** 2
    return np.sqrt(np.sum(diff2) / source_grad.shape[0])


def rmsd(history_file, reference_grad=None, stride=1.0):
    list_of_data = read_history_file(history_file, detect_rows(history_file))
    if reference_grad is not None:
        ref_data = pd.read_csv(reference_grad, delimiter=r'\s+', comment='#', header=None)
    else:
        ref_data = list_of_data[-1].copy()
    result = list()
    for i, src in enumerate(list_of_data, start=1):
        time = i * stride
        r = compute_rmsd(src, ref_data)
        result.append([time, r])
    return pd.DataFrame(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute gradient RMSD from history file and reference file')
    required_args = parser.add_argument_group('required named arguments')
    required_args.add_argument('--history', help='history file', required=True)
    required_args.add_argument('--output', help='output file', required=True)
    parser.add_argument('--reference', help='reference file')
    parser.add_argument('--stride', type=float, default=1.0, help='time stride of recording history')
    args = parser.parse_args()
    rmsd_result = rmsd(args.history, args.reference, args.stride)
    with open(args.output, 'w') as f_output:
        f_output.write('# time(ns) RMSD\n')
        rmsd_result.to_string(f_output, float_format='%15.10e', index=False, header=False)
        f_output.write('\n')

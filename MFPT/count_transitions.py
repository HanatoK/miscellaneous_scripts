#!/usr/bin/env python3
import pandas as pd
import numpy as np
# import regex as re
import csv
from math import isclose


def direct_find(s: str, start, end):
    ids = []
    current_start = None
    current_end = None
    for i, ch in enumerate(s):
        if ch == start and current_start is None:
            current_start = i
        if ch == end and current_end is None and current_start is not None:
            current_end = i
        if current_start is not None and current_end is not None:
            ids.append((current_start, current_end+1))
            current_start = None
            current_end = None
    return ids


def label_to_char(x):
    if isclose(x, 1.0):
        return 'B'
    elif isclose(x, 0.0):
        return 'A'
    else:
        return 'M'


def main():
    traj = pd.read_csv('traj.csv.gz')
    data = traj[['step', 'label']].to_numpy()
    # map the labels to a string
    labels_num = data[:, 1]
    f = np.vectorize(label_to_char)
    labels_str = ''.join(list(f(labels_num)))
    print('Finished reading')
    print('Start computing k_AB')
    # find all MFP subtrajs using regex
    # pattern = re.compile('A[AM]*B')
    # transition_A_to_B = [m.span() for m in re.finditer(pattern, labels_str)]
    transition_A_to_B = direct_find(labels_str, 'A', 'B')
    print(f'Number of transition_A_to_B: {len(transition_A_to_B)}')
    time_AB = []
    with open('transition_A_to_B.csv', 'w', newline='') as f_output:
        fieldnames = ['step_A', 'step_B', 'label_A', 'label_B']
        writer = csv.writer(f_output)
        writer.writerow(fieldnames)
        for i, j in transition_A_to_B:
            time_AB.append(data[j-1][0] - data[i][0])
            writer.writerow([data[i][0], data[j-1][0], data[i][1], data[j-1][1]])
    time_AB = np.array(time_AB) * 0.5 / 1e3
    print(time_AB)
    avg_time_AB = np.average(time_AB)
    stddev_time_AB = np.std(time_AB)
    rate_AB = 1.0/avg_time_AB
    stddev_rate_AB = np.sqrt(np.square(stddev_time_AB / avg_time_AB)) * rate_AB
    print(f'Mean first passage time from A to B is {avg_time_AB} ± {stddev_time_AB} ps, rate k_ab is {rate_AB:.3e} ps^-1 ± {stddev_rate_AB:.3e} ps^-1')
    # pattern = re.compile('B[BM]*A')
    print('Start computing k_BA')
    # transition_B_to_A = [m.span() for m in re.finditer(pattern, labels_str)]
    transition_B_to_A = direct_find(labels_str, 'B', 'A')
    print(f'Number of transition_B_to_A: {len(transition_B_to_A)}')
    time_BA = []
    with open('transition_B_to_A.csv', 'w', newline='') as f_output:
        fieldnames = ['step_B', 'step_A', 'label_B', 'label_A']
        writer = csv.writer(f_output)
        writer.writerow(fieldnames)
        for i, j in transition_B_to_A:
            time_BA.append(data[j-1][0] - data[i][0])
            writer.writerow([data[i][0], data[j-1][0], data[i][1], data[j-1][1]])
    time_BA = np.array(time_BA) * 0.5 / 1e3
    avg_time_BA = np.average(time_BA)
    stddev_time_BA = np.std(time_BA)
    rate_BA = 1.0/avg_time_BA
    stddev_rate_BA = np.sqrt(np.square(stddev_time_BA / avg_time_BA)) * rate_BA
    print(f'Mean first passage time from B to A is {avg_time_BA} ± {stddev_time_BA} ps, rate k_ba is {rate_BA:.3e} ps^-1 ± {stddev_rate_BA:.3e} ps^-1')


if __name__ == '__main__':
    main()

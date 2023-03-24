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
    with open('transition_A_to_B.csv', 'w', newline='') as f_output:
        time_avg = 0
        fieldnames = ['step_A', 'step_B', 'label_A', 'label_B']
        writer = csv.writer(f_output)
        writer.writerow(fieldnames)
        for i, j in transition_A_to_B:
            time_avg += data[j-1][0] - data[i][0]
            writer.writerow([data[i][0], data[j-1][0], data[i][1], data[j-1][1]])
        time_avg *= 0.5 / 1e3 / len(transition_A_to_B)
        print(f'Mean first passage time from A to B is {time_avg} ps, rate k_ab is {1.0/time_avg:.3e} ps^-1')
    # pattern = re.compile('B[BM]*A')
    print('Start computing k_BA')
    # transition_B_to_A = [m.span() for m in re.finditer(pattern, labels_str)]
    transition_B_to_A = direct_find(labels_str, 'B', 'A')
    print(f'Number of transition_B_to_A: {len(transition_B_to_A)}')
    with open('transition_B_to_A.csv', 'w', newline='') as f_output:
        time_avg = 0
        fieldnames = ['step_B', 'step_A', 'label_B', 'label_A']
        writer = csv.writer(f_output)
        writer.writerow(fieldnames)
        for i, j in transition_B_to_A:
            time_avg += data[j-1][0] - data[i][0]
            writer.writerow([data[i][0], data[j-1][0], data[i][1], data[j-1][1]])
        time_avg *= 0.5 / 1e3 / len(transition_B_to_A)
        print(f'Mean first passage time from B to A is {time_avg} ps, rate k_ba is {1.0/time_avg:.3e} ps^-1')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import csv
from math import isclose


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
    # find all MFP subtrajs using regex
    transition_A_to_B = [m.span() for m in re.finditer('A[A|M]*B', labels_str)]
    transition_B_to_A = [m.span() for m in re.finditer('B[B|M]*A', labels_str)]
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
    with open('transition_B_to_A.csv', 'w', newline='') as f_output:
        time_avg = 0
        fieldnames = ['step_B', 'step_A', 'label_B', 'label_A']
        writer = csv.writer(f_output)
        writer.writerow(fieldnames)
        for i, j in transition_B_to_A:
            time_avg += data[j-1][0] - data[i][0]
            writer.writerow([data[i][0], data[j-1][0], data[i][1], data[j-1][1]])
        time_avg *= 0.5 / 1e3 / len(transition_A_to_B)
        print(f'Mean first passage time from B to A is {time_avg} ps, rate k_ba is {1.0/time_avg:.3e} ps^-1')


if __name__ == '__main__':
    main()

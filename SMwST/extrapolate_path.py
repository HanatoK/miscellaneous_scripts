#!/usr/bin/env python3
import numpy as np

def extrapolatePath(pathFile):
    data = np.genfromtxt(pathFile, unpack=True)
    result = []
    tmp_list = []
    for cv in data:
        tmp_list.clear()
        # extrapolate the start point
        x = cv[0] - (cv[1] - cv[0])
        tmp_list.append(x)
        tmp_list.extend(cv.tolist())
        # extrapolate the end point
        x = cv[-1] - (cv[-2] - cv[-1])
        tmp_list.append(x)
        result.append(np.array(tmp_list))
    return np.array(result)

file_list = ['path_10copies_10ps.txt.origin', 'path_20copies_10ps.txt.origin',
             'path_20copies_50ps.txt.origin', 'path_20copies_100ps.txt.origin',
             'path_50copies_10ps.txt.origin']

for path_file in file_list:
    data = extrapolatePath(path_file)
    new_filename = path_file[:len(path_file)-len('.origin')]
    np.savetxt(new_filename, data.transpose(), fmt='%10.5f')

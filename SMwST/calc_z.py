#!/usr/bin/env python3
import numpy as np

def get_calc_func(pathFile):
    data = np.genfromtxt(pathFile, unpack=True)
    num_images = len(data[0])
    dcv = np.diff(data, axis=1)
    distances = np.sqrt(np.sum(dcv * dcv, axis=0))
    p_lambda = 1.0 / np.mean(distances * distances)
    print(f'{pathFile}: lambda = {p_lambda:12.7f} ; average_distance = {np.mean(distances):12.7f} ; variance = {np.var(distances):12.7f}')
    def calc_s(point):
        index_list = np.arange(0, num_images, 1)
        dist = -p_lambda * (point[:,None] - data) * (point[:,None] - data)
        numerator = np.sum(index_list * np.exp(np.sum(dist, axis=0)))
        denominator = np.sum(np.exp(np.sum(dist, axis=0)))
        factor = 1.0 / (num_images - 1)
        s = numerator / denominator * factor
        return s
    def calc_z(point):
        dist = -p_lambda * (point[:,None] - data) * (point[:,None] - data)
        dist = np.sum(np.exp(np.sum(dist, axis=0)))
        z = -1.0 / p_lambda * np.log(dist)
        return z
    return (calc_s, calc_z)

input_files = ['path_20copies_10ps.txt', 'path_10copies_10ps.txt', 
               'path_20copies_50ps.txt', 'path_20copies_100ps.txt',
               'path_20copies_100ps_2.txt']
calc_func_list = map(get_calc_func, input_files)

for pathFile, func in zip(input_files, calc_func_list):
    data = np.genfromtxt(pathFile, unpack=True).transpose()
    print(f'Path file: {pathFile}')
    for index, point in enumerate(data):
        print(f'Index = {index:3d} ; s = {func[0](point):12.7f} ; z = {func[1](point):12.7f}')

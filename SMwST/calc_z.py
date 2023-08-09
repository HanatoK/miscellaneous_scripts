#!/usr/bin/env python3
import numpy as np
from scipy.special import logsumexp, softmax


def numerical_grad(func, x, eps=0.001):
    num_grad = []
    for i in range(len(x)):
        prev_x = x.copy()
        next_x = x.copy()
        prev_x[i] -= eps
        next_x[i] += eps
        num_grad.append((func(next_x) - func(prev_x)) / (2.0 * eps))
    return np.array(num_grad)


def get_calc_func(pathFile):
    data = np.genfromtxt(pathFile, unpack=True)
    num_images = len(data[0])
    dcv = np.diff(data, axis=1)
    distances = np.sqrt(np.sum(dcv * dcv, axis=0))
    p_lambda = 1.0 / np.mean(distances * distances)
    print(f'{pathFile}: lambda = {p_lambda:12.7f} ; average_distance = {np.mean(distances):12.7f} ; variance = {np.var(distances):12.7f}')

    def calc_s(point):
        index_list = np.arange(0, num_images, 1)
        dist = -p_lambda * (point[:, None] - data) * (point[:, None] - data)
        s = 1.0 / (num_images - 1) * np.sum(index_list * softmax(np.sum(dist, axis=0)))
        return s

    def calc_s_derivative(point):
        B = calc_s(point) * (num_images - 1)
        index_list = np.arange(0, num_images, 1)
        diff = (point[:, None] - data)
        dist = -p_lambda * diff * diff
        Ai = softmax(np.sum(dist, axis=0))
        dsdx = -2.0 * p_lambda * diff / (num_images - 1) * (index_list - B) * Ai
        return dsdx

    def calc_z(point):
        dist = -p_lambda * (point[:, None] - data) * (point[:, None] - data)
        z = -1.0 / p_lambda * logsumexp(np.sum(dist, axis=0))
        return z

    def calc_z_derivative(point):
        diff = (point[:, None] - data)
        dist = -p_lambda * diff * diff
        Ai = softmax(np.sum(dist, axis=0))
        return 2.0 * diff * Ai

    return (calc_s, calc_z, calc_s_derivative, calc_z_derivative)


input_files = ['remove_loops/example/path_copy_19_iter_003.txt']
calc_func_list = map(get_calc_func, input_files)

for pathFile, func in zip(input_files, calc_func_list):
    data = np.genfromtxt(pathFile, unpack=True).transpose()
    print(f'Path file: {pathFile}')
    print('Calculate PCV s and z:')
    for index, point in enumerate(data):
        print(f'Index = {index:3d} ; s = {func[0](point):12.7f} ; z = {func[1](point):12.7f}')
    print('Verify derivatives of PCV s and z:')
    for index, point in enumerate(data):
        dsdx = np.sum(func[2](point), axis=1)
        dzdx = np.sum(func[3](point), axis=1)
        eps = 0.1
        for i in range(4):
            dsdx_num = numerical_grad(func[0], point, eps)
            dzdx_num = numerical_grad(func[1], point, eps)
            ds_err = np.sqrt(np.sum(np.mean(np.square(dsdx - dsdx_num))))
            dz_err = np.sqrt(np.sum(np.mean(np.square(dzdx - dzdx_num))))
            print(f'Index = {index:3d} ; epsilon = {eps:12.7e} ; ds_err = {ds_err:12.7e} ; dz_err = {dz_err:12.7e}')
            eps /= 10.0

#!/usr/bin/env python3
import numpy as np

def circular_mean(angular_data, weights=None):
    if weights is None:
        weights = np.ones(angular_data.shape[0])
    complex_data = np.exp(1.0j * angular_data)
    if len(complex_data.shape) == 1:
        complex_mean = np.sum(complex_data * weights, axis=0)
    elif len(weights.shape) > 1 or (weights.shape[0] != angular_data.shape[0]):
        raise RuntimeError(f'weights has shape {weights.shape}, but shape ({angular_data.shape[0]}) is required')
    else:
        complex_mean = np.sum(complex_data * weights[:, None], axis=0)
    return np.angle(complex_mean)

def circular_variance(angular_data, weights=None):
    if weights is None:
        weights = np.ones(angular_data.shape[0])
    complex_data = np.exp(1.0j * angular_data)
    if len(complex_data.shape) == 1:
        complex_mean = np.sum(complex_data * weights, axis=0) / np.sum(weights)
    elif len(weights.shape) > 1 or (weights.shape[0] != angular_data.shape[0]):
        raise RuntimeError(f'weights has shape {weights.shape}, but shape ({angular_data.shape[0]}) is required')
    else:
        complex_mean = np.sum(complex_data * weights[:, None], axis=0) / np.sum(weights)
    R_mean = np.absolute(complex_mean)
    return 1.0 - R_mean

def circular_stddev(angular_data, weights=None):
    complex_data = np.exp(1.0j * angular_data)
    if weights is None:
        weights = np.ones(angular_data.shape[0])
    elif len(weights.shape) > 1 or (weights.shape[0] != angular_data.shape[0]):
        raise RuntimeError(f'weights has shape {weights.shape}, but shape ({angular_data.shape[0]}) is required')
    if len(complex_data.shape) == 1:
        complex_mean = np.sum(complex_data * weights, axis=0) / np.sum(weights)
    else:
        complex_mean = np.sum(complex_data * weights[:, None], axis=0) / np.sum(weights)
    R_mean = np.absolute(complex_mean)
    return np.sqrt(-2.0 * np.log(R_mean))

def circular_correlation(theta_true, theta_pred, weights=None):
    if weights is None:
        weights = np.ones(theta_pred.shape[0])
    theta_pred_mean = circular_mean(theta_pred, weights)
    theta_true_mean = circular_mean(theta_true, weights)
    if len(theta_pred.shape) > 1:
        weights = weights[:, None]
    theta_pred_sin = np.sin(theta_pred - theta_pred_mean)
    theta_true_sin = np.sin(theta_true - theta_true_mean)
    tmp1 = np.sum(theta_pred_sin * theta_true_sin * weights, axis=0) / np.sum(weights)
    tmp2 = np.sum(theta_true_sin * theta_true_sin * weights, axis=0) / np.sum(weights)
    tmp3 = np.sum(theta_pred_sin * theta_pred_sin * weights, axis=0) / np.sum(weights)
    return tmp1 / np.sqrt(tmp2 * tmp3)

if __name__ == '__main__':
    # some tests
    def test_1():
        print(f'Test 1: mean, variance and standard deviation')
        angles = np.array([-165.0, 175.0, 60.0, -90.0])
        weights = np.array([1.0, 2.0, 1.0, 1.0])
        mean_1 = circular_mean(np.radians(angles), weights=weights)
        var_1 = circular_variance(np.radians(angles), weights=weights)
        stddev_1 = circular_stddev(np.radians(angles), weights=weights)
        print(f'angles = {angles} ; weights = {weights}\nmean = {mean_1} ; var = {var_1} ; stddev = {stddev_1}')
        angles = np.array([-165.0, 175.0, 60.0, -90.0, 175.0])
        weights = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        mean_2 = circular_mean(np.radians(angles), weights=weights)
        var_2 = circular_variance(np.radians(angles), weights=weights)
        stddev_2 = circular_stddev(np.radians(angles), weights=weights)
        print(f'angles = {angles} ; weights = {weights}\nmean = {mean_2} ; var = {var_2} ; stddev = {stddev_2}')
    def test_2():
        angles = np.array([-165.0, 175.0, 60.0, -90.0])
        weights = np.array([1.0, 2.0, 1.0, 1.0])
        angles_pred = np.array([-168.0, 179.0, 55.0, -93.0])
        corr = circular_correlation(theta_true=np.radians(angles),
                                    theta_pred=np.radians(angles_pred),
                                    weights=weights)
        print(f'Test 2: correlation')
        print(f'angles_true = {angles}\nangles_pred = {angles_pred}\ncorrelation coefficient = {corr}')
        angles = np.array([-165.0, 175.0, 60.0, -90.0, 175.0])
        weights = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        angles_pred = np.array([-168.0, 179.0, 55.0, -93.0, 179.0])
        corr = circular_correlation(theta_true=np.radians(angles),
                                    theta_pred=np.radians(angles_pred),
                                    weights=weights)
        print(f'angles_true = {angles}\nangles_pred = {angles_pred}\ncorrelation coefficient = {corr}')
    test_1()
    test_2()

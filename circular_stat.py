#!/usr/bin/env python3
import numpy as np

def circular_mean(angular_data, weights=None):
    if weights is None:
        weights = np.ones(angular_data.shape[0])
    complex_data = np.exp(1.0j * angular_data)
    if len(complex_data.shape) == 1:
        complex_mean = np.sum(complex_data * weights, axis=0)
    else:
        complex_mean = np.sum(complex_data * weights[:, None], axis=0)
    return np.angle(complex_mean)

def circular_variance(angular_data, weights=None):
    if weights is None:
        weights = np.ones(angular_data.shape[0])
    complex_data = np.exp(1.0j * angular_data)
    if len(complex_data.shape) == 1:
        complex_mean = np.sum(complex_data * weights, axis=0) / np.sum(weights)
    else:
        complex_mean = np.sum(complex_data * weights[:, None], axis=0) / np.sum(weights)
    R_mean = np.absolute(complex_mean)
    return 1.0 - R_mean

# equivalent to circstd(angular_data, weights=weights, method='circular') in astropy
def circular_stddev(angular_data, weights=None):
    complex_data = np.exp(1.0j * angular_data)
    if weights is None:
        weights = np.ones(angular_data.shape[0])
    if len(complex_data.shape) == 1:
        complex_mean = np.sum(complex_data * weights, axis=0) / np.sum(weights)
    else:
        complex_mean = np.sum(complex_data * weights[:, None], axis=0) / np.sum(weights)
    R_mean = np.absolute(complex_mean)
    return np.sqrt(-2.0 * np.log(R_mean))

if __name__ == '__main__':
    test_angles = np.array([-165.0, 175.0, 175.0])
    test_angles = np.radians(test_angles)
    print(circular_mean(test_angles))
    print(circular_variance(test_angles))
    print(circular_stddev(test_angles))
    test_angles = np.array([-165.0, 175.0])
    test_angles = np.radians(test_angles)
    test_weights = np.array([1.0,2])
    print(circular_mean(test_angles, test_weights))
    print(circular_variance(test_angles, test_weights))
    print(circular_stddev(test_angles, test_weights))

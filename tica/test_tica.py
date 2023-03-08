#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from scipy.linalg import eigh


def time_correlation_matrix(data, weights=None, tau=10):
    dim = np.shape(data)[1]
    length = np.shape(data)[0]
    if weights is None:
        weights = np.ones(length)
    weighted_data = data * weights[:, None]
    mean = np.sum(weighted_data, axis=0) / np.sum(weights)
    mat_C = np.zeros((dim, dim))
    mean_free_data = data - mean
    X_0 = mean_free_data[:length-tau] * weights[:length-tau, None]
    X_t = mean_free_data[tau:] * weights[tau:, None]
    mat_C = np.matmul(X_0.T, X_t) / np.sum(weights[tau:, None] * weights[:length-tau, None])
    mat_C = 0.5 * (mat_C + mat_C.T)
    return mat_C


def correlation_fftn(X):
    length = np.shape(X)[0]
    factors = np.arange(length, 0, -1)
    dim = np.shape(X)[1]
    mean = np.sum(X, axis=0) / length
    mean_free_data = X - mean
    mean_free_data = np.pad(mean_free_data, ((length-1, 0), (0, 0)),
                            mode='constant', constant_values=0)
    fft_X = np.fft.fft2(mean_free_data, axes=[0])
    fft_Y = np.fft.ifft2(mean_free_data, axes=[0])
    L = list()
    for i in range(0, dim):
        tmp = list()
        for j in range(0, dim):
            tmp.append(np.fft.fft(fft_X[:, i] * fft_Y[:, j])[:length] / factors)
        L.append(tmp)
    return np.transpose(np.array(L), axes=(2, 0, 1))


def tica(data, time_lag):
    cov_mat_t = time_correlation_matrix(data, tau=time_lag)
    cov_mat_0 = time_correlation_matrix(data, tau=0)
    return eigh(a=cov_mat_t, b=cov_mat_0)


def implied_time_scale(time_lag, eigenvalue):
    return -1.0 * time_lag / np.log(np.abs(eigenvalue))


def main():
    data = np.load('hmm-doublewell-2d-100k.npz')
    traj = data['trajectory']
    vals = np.ones(len(traj))
    # histogram data
    stats, x_edge, y_edge, binnumber = binned_statistic_2d(
        x=traj[:, 0], y=traj[:, 1], values=vals, statistic='sum',
        bins=(100, 50))
    # Histogram does not follow the Cartesian convention
    stats = np.transpose(stats)
    plt.imshow(stats, extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]])
    # tica
    time_lags = np.arange(1, 31)
    for time_lag in time_lags:
        eigvals, eigvecs = tica(traj, time_lag)
        its = implied_time_scale(time_lag=time_lag, eigenvalue=eigvals[-1])
        print(f'Time lag = {time_lag}; its = {its}')
        max_eigvec = eigvecs[:, -1] * 1.5
        plt.plot([0, max_eigvec[0]], [-1.0, max_eigvec[1]], linewidth=3)
    plt.show()


if __name__ == '__main__':
    main()


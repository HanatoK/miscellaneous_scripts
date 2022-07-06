#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


if __name__ == '__main__':
    with np.load('hmm-doublewell-2d-100k.npz') as fh:
        data = pd.DataFrame({'x': fh['trajectory'][:, 0],
                             'y': fh['trajectory'][:, 1]})
    cov_mat = time_correlation_matrix(data.to_numpy(), tau=10)
    print(cov_mat)
    eigvals, eigvecs = np.linalg.eigh(cov_mat)
    print(eigvecs)
    print(eigvals)
    plt.hist2d(data['x'], data['y'], bins=(50, 50), cmap=plt.cm.jet)
    max_eigvec = eigvecs[:, -1] * 3.0
    plt.plot([0, max_eigvec[0]], [0, max_eigvec[1]], linewidth=3)
    plt.show()

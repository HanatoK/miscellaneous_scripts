#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy import stats

class multivariate_linear_regression:

    def __init__(self, X, Y):
        import numpy as np
        # X and Y are 2D vectors
        # X.shape[0] is the number of points
        # X.shape[1] is the dimension
        # construct linear equation Ax=b
        # build A
        A = np.insert(X, 0, 1.0, axis=1)
        # helper function for solving Ax=b using projection
        def proj(p_A):
            return np.matmul(np.linalg.inv(np.matmul(p_A.T, p_A)),
                             p_A.T)
        self.factors = np.matmul(proj(A), Y)
        # print(self.factors)

    def get_K(self):
        return self.factors[1:]

    def get_b(self):
        return self.factors[0]

    def get_r2_score(self, X, y_true, sample_weight=None, mean_function=None, diff2_function=None):
        if sample_weight is None:
            sample_weight = np.ones(X.shape[0])
        if mean_function is None:
            def default_mean_function(data, sample_weight):
                return np.sum(data * np.transpose([sample_weight]), axis=0) / np.sum(sample_weight)
            mean_function = default_mean_function
        if diff2_function is None:
            def default_diff2(X1, X2, sample_weight):
                diff = X1 - X2
                return np.sum(np.square(diff) * np.transpose([sample_weight]), axis=0) / np.sum(sample_weight)
            diff2_function = default_diff2
        y_pred = self.evaluate(X)
        u = diff2_function(y_pred, y_true, sample_weight)
        y_true_mean = mean_function(y_true, sample_weight)
        v = diff2_function(y_true, y_true_mean, sample_weight)
        total_r2 = 1.0 - np.sum(u) / np.sum(v)
        component_r2 = 1.0 - u / v
        return total_r2, component_r2

    def evaluate(self, X):
        # evaluate a new dataset X
        A = np.insert(X, 0, 1.0, axis=1)
        return np.matmul(A, self.factors)

class multivariate_linear_regression_2:

    def __init__(self, X, Y):
        Y_T = Y.transpose()
        factors = list()
        for l in range(0, Y.shape[1]):
            factors.append(self._regress_1d(X, Y_T[l]))
        self.factors = np.transpose(factors)
        #print(self.factors)

    def _regress_1d(self, X, y):
        # y: n-elements vector
        # d: number of regressors
        n = X.shape[0]
        d = X.shape[1]
        # A: dxd matrix
        # b: d-elements column vector
        #A = np.zeros((d, d))
        b = np.zeros(d)
        X_T = X.transpose()
        K1 = np.full((n, n), -1.0)
        np.fill_diagonal(K1, n - 1.0)
        A = np.matmul(X_T, np.matmul(K1, X))
        for p in range(0, d):
            b[p] = n * np.sum(X_T[p] * y) - np.sum(y) * np.sum(X_T[p])
            #for q in range(p, d):
                #A[p][q] = n * np.sum(X_T[p] * X_T[q]) - np.sum(X_T[p]) * np.sum(X_T[q])
                #if q > p:
                    #A[q][p] = A[p][q]
        m = np.linalg.solve(A, b)
        m_0 = (np.sum(y) - np.sum(np.dot(X, m))) / n
        return np.insert(m, 0, m_0)

    def get_K(self):
        return self.factors[1:]

    def get_b(self):
        return self.factors[0]

# multivariate linear regression by gradient descent
class multivariate_linear_regression_sklearn:

    def __init__(self, X, Y):
        from sklearn import linear_model
        self.ols = linear_model.LinearRegression()
        self.model = self.ols.fit(X, Y)
        self.factors = np.insert(self.model.coef_.T, 0, self.model.intercept_, axis=0).copy()

    def get_K(self):
        return self.factors[1:]

    def get_b(self):
        return self.factors[0]

    def get_r2_score(self, X, y_true, sample_weight=None):
        y_pred = self.evaluate(X)
        from sklearn.metrics import r2_score
        total_r2 = r2_score(y_true, y_pred, sample_weight=sample_weight, multioutput='variance_weighted')
        component_r2 = r2_score(y_true, y_pred, sample_weight=sample_weight, multioutput='raw_values')
        return total_r2, component_r2

    def evaluate(self, X):
        # evaluate a new dataset X
        A = np.insert(X, 0, 1.0, axis=1)
        return np.matmul(A, self.factors)

def dihedral_mean(data, sample_weight):
    sin_data = np.sin(data)
    cos_data = np.cos(data)
    tmp1 = np.sum(np.arctan2(sin_data, cos_data) * np.transpose([sample_weight]), axis=0) / np.sum(sample_weight)
    return tmp1

def dihedral_diff2(X1, X2, sample_weight):
    # period is from -pi to pi
    def wrap(x, lower, upper):
        import numpy as np
        period = upper - lower
        if x >= lower and x <= upper:
            return x
        if x < lower:
            dist = lower - x
            num_period = np.floor(dist / period)
            tmp = np.abs(dist / period - np.round(dist / period))
            if np.isclose(tmp, 0):
                x += num_period * period
            else:
                x += (num_period + 1) * period
        if x > upper:
            dist = x - upper
            num_period = np.floor(dist / period)
            tmp = np.abs(dist / period - np.round(dist / period))
            if np.isclose(tmp, 0):
                x -= num_period * period
            else:
                x -= (num_period + 1) * period
        return x

    def periodic_dist(x, ref, lower, upper):
        import numpy as np
        vec_wrap = np.vectorize(wrap)
        x = vec_wrap(x, lower, upper)
        ref = vec_wrap(ref, lower, upper)
        dist = x - ref
        period = upper - lower
        length = 0.5 * period
        return np.where(np.abs(dist)>length, np.where(dist<0, dist+period, dist-period), dist)
        # if np.abs(dist) > length:
        #     if dist < 0:
        #         dist += period
        #     else:
        #         dist -= period
        # return dist

    diff = periodic_dist(X1, X2, -np.pi, np.pi)
    return np.sum(np.square(diff) * np.transpose([sample_weight]), axis=0) / np.sum(sample_weight)

if __name__ == '__main__':
    import numpy as np
    X = np.array([[-1.0, 1.5],
              [-6.2, -6.5],
              [2.3, 1.9],
              [-1.0, 1.0],
              [7.6, 1.8]])
    Y = np.array([[2.0, -4.3, 0.5, -1.4, -6.9]]).transpose()
    lr_model = multivariate_linear_regression(X, Y)
    print(X)
    print(Y)
    print(lr_model.factors)
    print(lr_model.evaluate(X))
    print(lr_model.get_r2_score(X, Y))
    lr_model = multivariate_linear_regression_2(X, Y)
    print(f'Using maximum likelihood:\n{lr_model.factors}')
    lr_model = multivariate_linear_regression_sklearn(X, Y)
    print(X)
    print(Y)
    print(lr_model.factors)
    print(lr_model.evaluate(X))
    print(lr_model.get_r2_score(X, Y))
    X = np.array([[-2, 0, 2,2], [-2, 1, 2,1]]).transpose()
    Y = np.array([[1,2,4,0], [1,3,4,-1]]).transpose()
    lr_model = multivariate_linear_regression(X, Y)
    print(X)
    print(Y)
    print(lr_model.factors)
    print(lr_model.evaluate(X))
    print(lr_model.get_r2_score(X, Y))
    lr_model = multivariate_linear_regression_sklearn(X, Y)
    print(X)
    print(Y)
    print(lr_model.factors)
    print(lr_model.evaluate(X))
    print(lr_model.get_r2_score(X, Y))
    lr_model = multivariate_linear_regression_2(X, Y)
    print(f'Using maximum likelihood:\n{lr_model.factors}')

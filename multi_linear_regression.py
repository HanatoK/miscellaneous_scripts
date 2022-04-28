#!/usr/bin/env python3

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

    def get_r2_score(self, X, y_true, sample_weight=None):
        import numpy as np
        if sample_weight is None:
            sample_weight = np.ones(X.shape[0])
        y_pred = self.evaluate(X)
        y_res = y_true - y_pred
        u = np.sum(np.square(y_res) * np.transpose([sample_weight]), axis=0) / np.sum(sample_weight)
        y_true_mean = np.sum(y_true * np.transpose([sample_weight]), axis=0) / np.sum(sample_weight)
        y_true_res = y_true - y_true_mean
        v = np.sum(np.square(y_true_res) * np.transpose([sample_weight]), axis=0) / np.sum(sample_weight)
        total_r2 = 1.0 - np.sum(u) / np.sum(v)
        component_r2 = 1.0 - u / v
        return total_r2, component_r2

    def evaluate(self, X):
        import numpy as np
        # evaluate a new dataset X
        A = np.insert(X, 0, 1.0, axis=1)
        return np.matmul(A, self.factors)

# multivariate linear regression by gradient descent
class multivariate_linear_regression_sklearn:

    def __init__(self, X, Y):
        from sklearn import linear_model
        ols = linear_model.LinearRegression()
        model = ols.fit(X, Y)
        self.factors = np.insert(model.coef_.T, 0, model.intercept_, axis=0).copy()

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
        import numpy as np
        # evaluate a new dataset X
        A = np.insert(X, 0, 1.0, axis=1)
        return np.matmul(A, self.factors)

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

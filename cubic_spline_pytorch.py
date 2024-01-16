#!/usr/bin/env python3
import torch


def cubic_spline_natural(X, Y):
    device = X.device
    N = len(X) - 1
    dX = X[1:] - X[:-1]
    dY = Y[1:] - Y[:-1]
    mat_rows = []
    vec_elems = [torch.zeros(1, device=device)]
    mat_rows.append(torch.cat([torch.tensor([1], device=device), torch.zeros(N, device=device)]))
    for i in range(1, N):
        zeros_before = torch.zeros(i-1, device=device)
        zeros_after = torch.zeros(N-2-(i-1), device=device)
        elem = torch.tensor([dX[i-1], 2.0 * (dX[i-1] + dX[i]), dX[i]], device=device)
        vec_elems.append(3.0 * (dY[i] / dX[i] - dY[i-1] / dX[i-1]))
        mat_rows.append(torch.cat([zeros_before, elem, zeros_after]))
    mat_rows.append(torch.cat([torch.zeros(N, device=device), torch.tensor([1], device=device)]))
    vec_elems.append(torch.zeros(1))
    mat = torch.stack(mat_rows)
    mat.to(device)
    vec = torch.as_tensor(vec_elems, device=device)
    tmp_C = torch.linalg.solve(mat, vec)
    m_A = Y[0:N]
    m_B = dY / dX - dX * (2.0 * tmp_C[0:N] + tmp_C[1:N+1]) / 3.0
    m_C = tmp_C[0:N]
    m_D = (tmp_C[1:N+1] - tmp_C[0:N]) / (3.0 * dX)

    def f(input_X):
        idx = torch.clamp(torch.searchsorted(X, input_X, right=True) - 1, min=0, max=N-1)
        idx.to(device)
        dx = input_X - X[idx]
        dx2 = dx * dx
        dx3 = dx2 * dx
        output_Y = m_A[idx] + m_B[idx] * dx + m_C[idx] * dx2 + m_D[idx] * dx3
        return output_Y

    def df(input_X):
        idx = torch.clamp(torch.searchsorted(X, input_X, right=True) - 1, min=0, max=N-1)
        idx.to(device)
        dx = input_X - X[idx]
        dx2 = dx * dx
        output_Y = m_B[idx] + 2.0 * m_C[idx] * dx + 3.0 * m_D[idx] * dx2
        return output_Y

    def df2(input_X):
        idx = torch.clamp(torch.searchsorted(X, input_X, right=True) - 1, min=0, max=N-1)
        idx.to(device)
        dx = input_X - X[idx]
        output_Y = 2.0 * m_C[idx] + 6.0 * m_D[idx] * dx
        return output_Y

    return f, df, df2


if __name__ == '__main__':

    def test_scipy():
        import numpy as np
        from scipy.interpolate import CubicSpline
        X = np.linspace(0, 5, 30)
        Y = np.cos(X)
        input_X = np.array([-0.1, 0.5, 1.5, 2.6, 3.7, 0.9, 5.2])
        fX = CubicSpline(X, Y, axis=0, bc_type='natural')
        output_Y = fX(input_X)
        print("Running test_scipy():")
        print(output_Y)

    def test_torch():
        X = torch.linspace(0, 5, 30, requires_grad=True)
        # X.requires_grad = True
        Y = torch.cos(X)
        # Y.requires_grad = True
        input_X = torch.tensor([-0.1, 0.5, 1.5, 2.6, 3.7, 0.9, 5.2], requires_grad=True)
        fX, df, df2 = cubic_spline_natural(X, Y)
        output_Y = fX(input_X)
        loss = torch.sum(output_Y)
        # loss.backward()
        grad_input_X = torch.autograd.grad(loss, input_X, create_graph=True)[0]
        grad2_input_X = torch.autograd.grad(torch.sum(grad_input_X), input_X)[0]
        print("Running test_torch():")
        print(output_Y)
        print(grad_input_X)
        print(df(input_X))
        print(grad2_input_X)
        print(df2(input_X))

    def test_torch2():
        from copy import deepcopy
        X = torch.linspace(0, 5, 30, requires_grad=True)
        Y = torch.cos(X)
        fX, df, df2 = cubic_spline_natural(X, Y)
        input_X = deepcopy(X)
        input_X.requires_grad = True
        output_Y = fX(input_X)
        loss = torch.sum(output_Y)
        loss.backward()
        print("Running test_torch2():")
        print(output_Y)
        print(X.grad)
        print(input_X.grad)
        print(df(input_X))

    def test_integral():
        import numpy as np
        p = torch.tensor(np.random.normal(0, 0.5, 30))
        X = torch.linspace(0, 5, 30, requires_grad=True)
        Y = torch.cos(X)
        fX, df, df2 = cubic_spline_natural(X, Y)

        def loss_func(X, inputs):
            new_X = X + p
            new_Y = torch.cos(new_X)
            new_fX, new_df, new_df2 = cubic_spline_natural(new_X, new_Y)
            return torch.sum((new_df(inputs) - df(inputs))**2)
        inputs = torch.linspace(0, 5, 10, requires_grad=True)
        L = loss_func(X, inputs)
        print("Running test_integral():")
        print(L)
        L.backward()
        print(X.grad)

        # def numerical_gradients(epsilon=0.001):
        #     results = torch.zeros_like(X, requires_grad=False)
        #     for i in range(len(X)):
        #         delta = torch.zeros_like(X, requires_grad=False)
        #         delta[i] = epsilon
        #         X_next = X

    test_scipy()
    test_torch()
    test_torch2()
    test_integral()

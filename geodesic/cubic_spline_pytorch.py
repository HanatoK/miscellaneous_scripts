#!/usr/bin/env python3
import torch
import numpy as np
from scipy.interpolate import CubicSpline


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
    A = []
    B = []
    C = []
    D = []
    for i in range(0, N):
        C.append(tmp_C[i])
        D.append((tmp_C[i+1] - tmp_C[i]) / (3.0 * dX[i]))
        B.append(dY[i] / dX[i] - dX[i] * (2.0 * tmp_C[i] + tmp_C[i+1]) / 3.0)
        A.append(Y[i])
    m_A = torch.as_tensor(A, device=device)
    m_B = torch.as_tensor(B, device=device)
    m_C = torch.as_tensor(C, device=device)
    m_D = torch.as_tensor(D, device=device)

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
        X = np.linspace(0, 5, 30)
        Y = np.cos(X)
        input_X = np.array([-0.1, 0.5, 1.5, 2.6, 3.7, 0.9, 5.2])
        fX = CubicSpline(X, Y, axis=0, bc_type='natural')
        output_Y = fX(input_X)
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
        input_X.reqiures_grad = True
        output_Y = fX(input_X)
        loss = torch.sum(output_Y)
        loss.backward()
        print(output_Y)
        print(X.grad)
        print(input_X.grad)
        print(df(input_X))

    test_scipy()
    test_torch()
    test_torch2()

#!/usr/bin/env python3
from cubic_spline_pytorch import cubic_spline_natural
import torch
import argparse
import pandas as pd


def reparametrize(positions, num_images=None, resolution=1000):
    # do interpolation
    X = torch.linspace(0, resolution, len(positions), device=positions.device)
    # iterate over columns
    fX_list = []
    df_list = []
    df2_list = []
    num_variables = positions.shape[1]
    if num_images is None:
        num_images = len(positions)
    for i in range(num_variables):
        X.requires_grad_(True)
        fX, _, _ = cubic_spline_natural(X, positions[:, i])
        def df(new_X):
            new_X.requires_grad_(True)
            output_Y = fX(new_X)
            loss = torch.sum(output_Y)
            return torch.autograd.grad(loss, new_X, create_graph=True)[0]
        def df2(new_X):
            new_X.requires_grad_(True)
            output_Y = df(new_X)
            loss = torch.sum(output_Y)
            return torch.autograd.grad(loss, new_X, create_graph=True)[0]
        fX_list.append(fX)
        df_list.append(df)
        df2_list.append(df2)
    new_X = torch.linspace(0, resolution, num_images * resolution, device=positions.device)
    new_Y = torch.stack([f(new_X) for f in fX_list]).T
    new_dY = torch.stack([f(new_X) for f in df_list]).T
    new_dY2 = torch.stack([f(new_X) for f in df2_list]).T
    # compute the total length
    adjacent_diff = torch.diff(new_Y, dim=0)
    all_lengths = torch.sqrt(torch.sum(adjacent_diff * adjacent_diff, dim=1))
    total_length = torch.sum(all_lengths)
    # search points on the spline
    adjacent_length = total_length / (num_images - 1)
    cumsum_l = torch.cumsum(all_lengths[:-1], dim=0)
    cumsum_adj_len = torch.arange(1, num_images - 1, device=positions.device) * adjacent_length
    idx = torch.searchsorted(cumsum_l, cumsum_adj_len) + 1
    # print(torch.unsqueeze(new_Y[0], 0))
    # print(new_Y[idx])
    # print(torch.unsqueeze(new_Y[-1], 0))
    result = torch.concatenate([torch.unsqueeze(new_Y[0], 0), new_Y[idx], torch.unsqueeze(new_Y[-1], 0)])
    derivative = torch.concatenate([torch.unsqueeze(new_dY[0], 0), new_dY[idx], torch.unsqueeze(new_dY[-1], 0)])
    derivative2 = torch.concatenate([torch.unsqueeze(new_dY2[0], 0), new_dY2[idx], torch.unsqueeze(new_dY2[-1], 0)])
    result.to(positions.device)
    derivative.to(positions.device)
    derivative2.to(positions.device)
    return result, derivative, derivative2


def remove_loops(positions):
    import numpy as np

    def image_distance(pos_x, pos_y):
        diff = np.asarray(pos_x) - np.asarray(pos_y)
        return np.sqrt(np.sum(diff * diff))
    results = positions.tolist()
    has_loop = False
    while True:
        has_loop = False
        for i in range(1, len(results) - 1):
            neighbor_dist = image_distance(results[i-1], results[i])
            j = i + 1
            while j < len(results):
                dist = image_distance(results[i-1], results[j])
                if dist < neighbor_dist:
                    has_loop = True
                    break
                j = j + 1
            if has_loop:
                print(f'Remove indexes from {i} to {j}')
                del results[i:j]
                break
        if not has_loop:
            break
    return torch.as_tensor(results, device=positions.device, dtype=positions.dtype)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='path position file')
    parser.add_argument('-n', '--num_images', type=int, help='number of desired images (default to the same images of the path)')
    parser.add_argument('-i', '--num_iterations', default=1, type=int, help='number of iterations')
    parser.add_argument('-o', '--output', help='output filename')
    args = parser.parse_args()
    positions = pd.read_csv(args.path, delimiter=r'\s+', comment='#', header=None).to_numpy()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    positions = torch.tensor(positions, dtype=torch.float, device=device)
    # Y = torch.tensor(positions, dtype=torch.float)
    Y = positions.clone().detach()
    Y.to(device)
    for i in range(0, args.num_iterations):
        Y, dY, dY2 = reparametrize(positions=Y, num_images=args.num_images)
    # write the output
    with open(args.output, 'w') as fOutput:
        tmp_Y = Y.cpu().detach().numpy()
        for point in tmp_Y:
            for elem in point:
                fOutput.write(f' {elem:15.10f}')
            fOutput.write('\n')
    # debug: show neighboring distances
    adjacent_diff = torch.diff(Y, dim=0)
    all_lengths = torch.sqrt(torch.sum(adjacent_diff * adjacent_diff, axis=1))
    for i, l in enumerate(all_lengths):
        print(f"Distance between image {i} and {i+1}: {l:12.5f}")
    # debug: show moved distances if possible
    if Y.shape == positions.shape:
        dist_Y_pos = Y - positions
        distances = torch.sqrt(torch.sum(dist_Y_pos * dist_Y_pos, dim=1))
        for i, l in enumerate(distances):
            print(f'Image {i} has moved: {l:12.5f}')
    # debug: show the derivatives
    print('Derivative from cubic splines:')
    print(dY)
    print('Second-order derivative from cubic splines:')
    print(dY2)


if __name__ == '__main__':
    main()

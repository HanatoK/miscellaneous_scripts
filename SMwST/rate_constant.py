#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import colorama
from colorama import Fore, Back, Style

# compute the Kramers rate constant from A to B
def compute_rate_constant(A_index, B_index, beta, x_star_index, X, pmf, diffusivity):
    # the first integral
    dx = X[1] - X[0]
    first_integrand = np.exp(-beta * pmf) * dx
    tmp_integral = np.cumsum(first_integrand)
    first_integral = tmp_integral[x_star_index]
    print('Integrand 1:')
    print(first_integrand)
    print('Integration 1:')
    print(tmp_integral)
    # the second integral
    second_integrand = np.exp(beta * pmf) / diffusivity * dx
    tmp_integral = np.cumsum(second_integrand)
    print('Integration 2:')
    print(tmp_integral)
    second_integral = tmp_integral[B_index] - tmp_integral[A_index]
    # debug
    print(f'Integral 1: {first_integral:12.7f} ; integral 2 {second_integral:12.7f}')
    return 1.0 / (first_integral * second_integral)

# find the nearest index of x
def find_index(x, X):
    tmp = np.abs(X - x)
    return np.argmin(tmp)

if __name__ == '__main__':
    colorama.init()
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pmf', help='PMF file')
    parser.add_argument('--input_diffusivity', help='diffusivity file')
    parser.add_argument('--beta', type=float, default=(1.0/(300 * 0.0019872041)), help='inverse temperature (1/kBT)')
    parser.add_argument('--cv_A', type=float, help='CV value of basin A')
    parser.add_argument('--cv_B', type=float, help='CV value of basin B')
    parser.add_argument('--cv_x_star', type=float, help='CV value of x*')
    args = parser.parse_args()
    print(Back.RED + Fore.WHITE + 'WARNING: please double check the calculation!' + Style.RESET_ALL)
    # read data
    data_pmf = pd.read_csv(args.input_pmf, delimiter='\s+', header=None, comment='#')
    data_diffusivity = pd.read_csv(args.input_diffusivity, delimiter='\s+', header=None, comment='#')
    # get the array of CV values
    X = data_pmf[0].to_numpy()
    # find the indexes
    A_index = find_index(args.cv_A, X)
    B_index = find_index(args.cv_B, X)
    x_star_index = find_index(args.cv_x_star, X)
    rate_constant = compute_rate_constant(
        A_index=A_index, B_index=B_index, beta=args.beta,
        x_star_index=x_star_index, X=X, pmf=data_pmf[1].to_numpy(),
        diffusivity=data_diffusivity[1].to_numpy())
    print(f'k = {rate_constant:15.7f}')

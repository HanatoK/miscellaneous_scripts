#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse

def get_minimum_location(filename):
    data = pd.read_csv(filename, header=None, delimiter='\s+', comment='#')
    last_col = data.columns[-1]
    idxmin = data[[last_col]].idxmin().iloc[0]
    min_location = np.array(data.iloc[idxmin])[0:-1]
    return min_location

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("pmf", help = "specify the PMF file")
    args = parser.parse_args()
    min_location = get_minimum_location(args.pmf)
    print(f'Minimum location: {min_location}')

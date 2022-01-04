#!/usr/bin/env python3
import pandas as pd

num_copies = 20
num_iterations = 100

if __name__ == '__main__':
    for i in range(0, num_iterations):
        data_sum = None
        for j in range(0, num_copies):
            filename = 'path_copy_' + str(j).zfill(2) + '_' + str(i).zfill(3) + '.dat'
            if data_sum is None:
                data_sum = pd.read_csv(filename, delimiter='\s+', header=None, comment='#')
            else:
                data_sum = data_sum + pd.read_csv(filename, delimiter='\s+', header=None, comment='#')
        data_avg = data_sum / num_copies
        output_filename = f'path_{i:03d}'
        with open(output_filename, 'w') as fOutput:
            for k in range(0, len(data_avg)):
                for l in range(0, len(data_avg.loc[k])):
                    fOutput.write(f' {data_avg.loc[k, l]:10.5f}')
                fOutput.write('\n')

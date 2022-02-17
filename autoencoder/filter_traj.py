#!/usr/bin/env python3
import pandas as pd

def write_dataframe(df:pd.DataFrame, filename:str):
    with open(filename, 'w') as f_output:
        df.to_string(f_output, float_format=' %15.10f', index=False, header=False)
        f_output.write('\n')

if __name__ == '__main__':
    traj_data = pd.read_csv('alad_round1.colvars.traj', delimiter='\s+', comment='#', header=None)
    # filter basin1
    basin1_data = traj_data.loc[(traj_data[1] > -0.62) & (traj_data[1] < -0.52)]
    write_dataframe(basin1_data, 'basin1.traj')
    # filter basin2
    basin2_data = traj_data.loc[(traj_data[1] > -0.05) & (traj_data[1] < 0.02)]
    write_dataframe(basin2_data, 'basin2.traj')
    # filter basin3
    basin3_data = traj_data.loc[(traj_data[1] > 0.64) & (traj_data[1] < 0.87)]
    write_dataframe(basin3_data, 'basin3.traj')


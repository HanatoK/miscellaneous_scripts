#!/usr/bin/env python3
from plot_colvars_traj import Colvars_traj
import numpy as np


def filter_traj(df, columns, point, radius, output_file):
    r = None
    for cv, x in zip(columns, point):
        if r is None:
            r = (df[cv] - x) * (df[cv] - x)
        else:
            r += (df[cv] - x) * (df[cv] - x)
    condition = (r < radius * radius)
    print(r)
    df.loc[condition].to_csv(output_file, index=False)
    return df[condition]


if __name__ == '__main__':
    cv_traj = Colvars_traj(['alad.colvars.traj']).as_pandas()
    columns = ['CV1', 'CV2']
    cv_traj[columns] = np.degrees(cv_traj[columns])
    #print(cv_traj)
    filter_traj(cv_traj, columns, [-1.14, 1.45], 0.2, 'basin1.csv')
    filter_traj(cv_traj, columns, [1.15, 1.55], 0.2, 'basin2.csv')
    filter_traj(cv_traj, columns, [2.19, -1.03], 0.2, 'basin3.csv')

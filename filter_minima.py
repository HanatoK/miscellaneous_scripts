#!/usr/bin/env python3
from plot_colvars_traj import Colvars_traj


def mark_2D_minimum(df, cv_list, points, margins, mark_index):
    if 'minima_mark' not in df:
        df['minima_mark'] = 0
    condition = None
    for cv, mar, p in zip(cv_list, margins, points):
        lower = p - mar * 0.5
        upper = p + mar * 0.5
        tmp = (df[cv] > lower) & (df[cv] < upper)
        if condition is None:
            condition = tmp
        else:
            condition = tmp & condition
    df.loc[condition, 'minima_mark'] = mark_index
    return df


if __name__ == '__main__':
    cv_traj = Colvars_traj(['alad.colvars.traj'])
    cv_traj_pd = cv_traj.as_pandas()
    minima_index = [1, 2, 3, 4, 5]
    cv_traj_pd = mark_2D_minimum(
        df=cv_traj_pd, cv_list=['CV1', 'CV2'],
        points=[-0.6, -0.35], margins=[0.1, 0.4],
        mark_index=minima_index[0])
    cv_traj_pd = mark_2D_minimum(
        df=cv_traj_pd, cv_list=['CV1', 'CV2'],
        points=[-0.58, -0.7], margins=[0.2, 0.4],
        mark_index=minima_index[0])
    cv_traj_pd = mark_2D_minimum(
        df=cv_traj_pd, cv_list=['CV1', 'CV2'],
        points=[0.45, -0.2], margins=[0.1, 0.4],
        mark_index=minima_index[1])
    cv_traj_pd = mark_2D_minimum(
        df=cv_traj_pd, cv_list=['CV1', 'CV2'],
        points=[0.65, -0.7], margins=[0.15, 0.25],
        mark_index=minima_index[2])
    cv_traj_pd = mark_2D_minimum(
        df=cv_traj_pd, cv_list=['CV1', 'CV2'],
        points=[0.65, 0.3], margins=[0.15, 0.2],
        mark_index=minima_index[3])
    cv_traj_pd = mark_2D_minimum(
        df=cv_traj_pd, cv_list=['CV1', 'CV2'],
        points=[-0.4, 0.25], margins=[0.15, 0.2],
        mark_index=minima_index[4])
    cv_traj_pd = mark_2D_minimum(
        df=cv_traj_pd, cv_list=['CV1', 'CV2'],
        points=[0.1, 0.38], margins=[0.2, 0.1],
        mark_index=minima_index[4])
    for i in minima_index:
        min_traj = cv_traj_pd[cv_traj_pd['minima_mark'] == i]
        min_traj.to_csv(f'minima_{i}.csv', index=False)

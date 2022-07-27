#!/usr/bin/env python3
from read_colvars_traj import ReadColvarsTraj
import csv


with open('eq_traj.csv', 'w') as f_output:
    writer = csv.writer(f_output)
    first_line = True
    with ReadColvarsTraj('eq_sin_cos.colvars.traj') as f_input:
        for line in f_input:
            if first_line:
                writer.writerow(line.keys())
                first_line = False
            else:
                writer.writerow([line.get(key) for key in line.keys()])

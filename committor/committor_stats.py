#!/usr/bin/env python3
import glob

def determineAB(traj_file):
    with open(traj_file, 'r') as fInput:
        for line in fInput:
            line = line.strip()
            if line.startswith('#'):
                continue
            fields = line.split()
            cv_s = float(fields[1])
            if cv_s < 0.36:
                return 'A'
            elif cv_s > 0.61:
                return 'B'
            else:
                continue
    return 'ERROR'


if __name__ == '__main__':
    structure_dirlist = sorted(glob.glob('gpath_[0-9][0-9][0-9][0-9]'))
    for structure in structure_dirlist:
        output_dir = structure + '/output_committor_gpath/'
        traj_filelist = sorted(glob.glob(output_dir + 'committor_gpath_out_*.colvars.traj'))
        num_A = 0.0
        num_B = 0.0
        for traj_file in traj_filelist:
            result = determineAB(traj_file)
            if result == 'A':
                num_A = num_A + 1.0
            elif result == 'B':
                num_B = num_B + 1.0
            else:
                raise RuntimeError('Cannot determine whether A or B.')
        total = num_A + num_B
        prob_A = num_A / total
        print(f'{structure} P_A = {prob_A:12.7f}')

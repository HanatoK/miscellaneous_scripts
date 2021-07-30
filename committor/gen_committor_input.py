#!/usr/bin/env python3
import string
import glob
import os
import shutil

if __name__ == '__main__':
    pdb_filenames = sorted(glob.glob('gpath_structures/gpath*.pdb'))
    num_simulations = 500
    for i, pdb_filename in enumerate(pdb_filenames):
        simulation_dirname = f'gpath_{i:04d}'
        if not os.path.exists(simulation_dirname):
            os.makedirs(simulation_dirname)
            output_dirname = simulation_dirname + '/output_committor_gpath'
            if not os.path.exists(output_dirname):
                os.makedirs(output_dirname)
        for j in range(num_simulations):
            index_string = f'{i:04d}_{j:03d}'
            namd_config_filename = f'{simulation_dirname}/committor_gpath_namd_{index_string}.namd'
            with open('committor_gpath_namd.template', 'r') as fInput:
                template_string = string.Template(fInput.read())
                config = template_string.safe_substitute(
                    pdb_filename=f'../{pdb_filename}',
                    output_filename=f'committor_gpath_out_{index_string}')
                with open(namd_config_filename, 'w') as fOutput:
                    fOutput.write(config)
                shutil.copyfile('committor_gpath_colvars.in', f'{simulation_dirname}/committor_gpath_colvars.in')
                shutil.copyfile('path_output.txt', f'{simulation_dirname}/path_output.txt')

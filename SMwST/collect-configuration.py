#!/usr/bin/env python3
from glob import glob
import argparse
import os
from shutil import copyfile

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('n', help='number of images')
args = parser.parse_args()
num_images = int(args.n)

# get the number of trajectories per swarm
output_basedir = 'output'
output_dirs = glob(f'{output_basedir}/*')
num_swarms = len(output_dirs) // num_images
print(f'Number of images : {num_images}')
print(f'Number of trajectories per image : {num_swarms}')
output_dirs.sort()
result_dirs = [output_dirs[(start*num_swarms):(start*num_swarms+num_swarms)] for start in range(num_images)]

# get the number of iterations
num_iterations = len(glob(result_dirs[0][0] + '/*.image*.iter*.coor'))
digits_iterations = len(str(num_iterations-1))
print(f'Number of iterations: {num_iterations}')

# make directories and copy images
base_directory = os.getcwd() + '/string_results'
if not os.path.exists(base_directory):
    os.makedirs(base_directory)
digits_images = len(str(num_images))
digits_copies = len(str(num_swarms))
for i_image in range(0, num_images):
    image_dirname = base_directory + '/image_' + str(i_image).zfill(digits_images)
    if not os.path.exists(image_dirname):
        os.makedirs(image_dirname)
    for j_copy in range(0, num_swarms):
        copy_dirname = image_dirname + '/copy_' + str(j_copy).zfill(digits_copies)
        if not os.path.exists(copy_dirname):
            os.makedirs(copy_dirname)
        result_dirname = result_dirs[i_image][j_copy]
        iteration_images = glob(result_dirname + '/*.image*.iter*.coor')
        iteration_images.sort()
        for i_iteration, image_name in enumerate(iteration_images):
            dest_name = copy_dirname + '/' + str(i_iteration).zfill(digits_iterations) + '.coor'
            print(f'Copying file from {image_name} to {dest_name}')
            copyfile(image_name, dest_name)

# collect the final path into a separate directory
print('Collect the images of the final path:')
final_path_directory = os.getcwd() + '/final_path'
if not os.path.exists(final_path_directory):
    os.makedirs(final_path_directory)
for j_copy in range(0, num_swarms):
    dest_dirname = final_path_directory + '/copy_' + str(j_copy).zfill(digits_copies)
    if not os.path.exists(dest_dirname):
        os.makedirs(dest_dirname)
    for i_image in range(0, num_images):
        image_dirname = base_directory + '/image_' + str(i_image).zfill(digits_images)
        copy_dirname = image_dirname + '/copy_' + str(j_copy).zfill(digits_copies)
        last_iteration = num_iterations - 1
        src_name = copy_dirname + '/' + str(last_iteration).zfill(digits_iterations) + '.coor'
        dest_name = dest_dirname + '/' + str(i_image).zfill(digits_images) + '.coor'
        print(f'Copying file from {src_name} to {dest_name}')
        copyfile(src_name, dest_name)

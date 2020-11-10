#!/usr/bin/env python3
from glob import glob
import argparse
import os
from shutil import copyfile

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('n', help = 'number of images')
args = parser.parse_args()
num_images = int(args.n)

# get the number of trajectories per swarm
output_dirs = glob('output/*')
num_swarms = int(len(output_dirs) / num_images)
print(f'Number of images : {num_images}')
print(f'Number of trajectories per image : {num_swarms}')
output_dirs.sort()
result_dirs = output_dirs[::num_swarms]

# get the number of iterations
num_iterations = len(glob(result_dirs[0] + '/*.image*.iter*.coor'))
digits_iterations = len(str(num_iterations))
print(f'Number of iterations: {num_iterations}')

# make directories and copy images
base_directory = os.getcwd() + '/string_results'
if not os.path.exists(base_directory):
    os.makedirs(base_directory)
i_image = 0
i_iteration = 0
digits_images = len(str(num_images))
for result_dirname in result_dirs:
    image_dirname = base_directory + '/' + str(i_image).zfill(digits_images)
    if not os.path.exists(image_dirname):
        os.makedirs(image_dirname)
    iteration_images = glob(result_dirname + '/*.image*.iter*.coor')
    iteration_images.sort()
    for image_name in iteration_images:
        dest_name = image_dirname + '/' + str(i_iteration).zfill(digits_iterations) + '.coor'
        print(f'Copying file from {image_name} to {dest_name}')
        copyfile(image_name, dest_name)
        i_iteration = i_iteration + 1
    i_image = i_image + 1
    i_iteration = 0

# collect the final path into a separate directory
print('Collect the images of the final path:')
final_path_directory = os.getcwd() + '/final_path'
if not os.path.exists(final_path_directory):
    os.makedirs(final_path_directory)
for i_image in range(0, num_images):
    image_dirname = base_directory + '/' + str(i_image).zfill(digits_images)
    last_iteration = num_iterations - 1
    last_image_name = image_dirname + '/' + str(last_iteration).zfill(digits_iterations) + '.coor'
    dest_name = final_path_directory + '/' + str(i_image).zfill(digits_images) + '.coor'
    print(f'Copying file from {last_image_name} to {dest_name}')
    copyfile(last_image_name, dest_name)

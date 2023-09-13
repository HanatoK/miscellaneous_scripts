#!/usr/bin/env python3
from glob import glob
import argparse
import os
from shutil import copyfile

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('n', help='number of images')
parser.add_argument('--output_dir', type=str, default='output/', help='output base directory')
args = parser.parse_args()
num_images = int(args.n)

output_basedir = args.output_dir
output_dirs = sorted(glob(f'{output_basedir}/*'))
print(f'Number of images : {num_images}')

# collect results of all iterations
# make directories and copy images
base_directory = os.getcwd() + '/string_results'
print(f'Collect the images of all iterations to {base_directory}')
if not os.path.exists(base_directory):
    os.makedirs(base_directory)
digits_images = len(str(num_images))
for i_image, i_res in zip(range(0, num_images), output_dirs):
    image_dirname = base_directory + '/image_' + str(i_image).zfill(digits_images)
    iteration_images = sorted(glob(i_res + '/*.image*.iter*.coor'))
    num_iterations = len(iteration_images)
    digits_iterations = len(str(num_iterations))
    if not os.path.exists(image_dirname):
        os.makedirs(image_dirname)
    for i_iteration, image_name in enumerate(iteration_images):
        dest_name = image_dirname + '/' + str(i_iteration).zfill(digits_iterations) + '.coor'
        print(f'Copying file from {image_name} to {dest_name}')
        copyfile(image_name, dest_name)

# collect the final path into a separate directory
final_path_directory = os.getcwd() + '/final_path'
print(f'Collect the images of the final path to {final_path_directory}')
if not os.path.exists(final_path_directory):
    os.makedirs(final_path_directory)
for i_image in range(0, num_images):
    image_dirname = base_directory + '/image_' + str(i_image).zfill(digits_images)
    iteration_images = sorted(glob(image_dirname + '/*.coor'))
    num_iterations = len(iteration_images)
    last_iteration = num_iterations - 1
    src_name = image_dirname + '/' + str(last_iteration).zfill(digits_iterations) + '.coor'
    dest_name = final_path_directory + '/' + str(i_image).zfill(digits_images) + '.coor'
    print(f'Copying file from {src_name} to {dest_name}')
    copyfile(src_name, dest_name)

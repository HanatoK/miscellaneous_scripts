#!/usr/bin/env python3
import numpy as np
import argparse


def get_path(start, end, num_images):
    num_middle_images = num_images[1]
    middle_path = np.linspace(start, end, num_middle_images)
    diff = np.diff(middle_path, axis=0)[0]
    prepend_path = []
    for i in range(num_images[0]):
        prepend_path.append(start - diff * (i + 1))
    append_path = []
    for i in range(num_images[2]):
        append_path.append(end + diff * (i + 1))
    path_list = []
    if prepend_path:
        path_list.append(np.array(list(reversed(prepend_path))))
    path_list.append(middle_path)
    if append_path:
        path_list.append(append_path)
    path = np.concatenate(path_list)
    return path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', nargs='+', type=float, help='start point')
    parser.add_argument('--end', nargs='+', type=float, help='end point')
    parser.add_argument(
        '--num_images', nargs=3, type=int,
        help='specify the number of images before the start point,'
        'between the start and end points, and after the end point')
    parser.add_argument('--output', type=str, help='output file')
    args = parser.parse_args()
    start_point = args.start
    end_point = args.end
    num_images = args.num_images
    path = get_path(np.array(start_point), np.array(end_point), num_images)
    if args.output:
        np.savetxt(args.output, path)
    print(f'Total number of images: {len(path)}')
    print(path)
    # print(np.diff(path, axis=0))

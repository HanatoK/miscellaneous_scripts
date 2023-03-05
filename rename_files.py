#!/usr/bin/env python3
import os
import shutil
from rich import print
from glob import glob
from natsort import natsorted


def rename_subtitles(video_files, subtitle_files):
    for video, subtitle in zip(video_files, subtitle_files):
        name_prefix = os.path.splitext(video)[0]
        subtitle_suffix = os.path.splitext(subtitle)[-1]
        new_name = name_prefix + subtitle_suffix
        print(f'Copy {subtitle} to {new_name}')
        shutil.copyfile(subtitle, new_name)


def main():
    video_files = natsorted(glob('* [0-9][0-9] (BDrip*.mkv'))
    subtitle_files = natsorted(glob('sc/*.ass'))
    print(video_files)
    print(subtitle_files)
    rename_subtitles(video_files, subtitle_files)


if __name__ == '__main__':
    main()

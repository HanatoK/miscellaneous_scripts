#!/bin/sh
dirlist=(gpath_[0-9][0-9][0-9][0-9])
for gpath_dir in ${dirlist[@]};
do
  echo "Compress $gpath_dir"
  tar cJf $gpath_dir.tar.xz $gpath_dir
done

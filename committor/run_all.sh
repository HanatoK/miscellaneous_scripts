#!/bin/bash
export LD_LIBRARY_PATH=/usr/local/lib64:/home/lia/hong/haochuan/lib64:$LD_LIBRARY_PATH
# namd_executable=/home/PrismXD/namd/cpu_build/Linux-x86_64-g++/namd2
namd_executable=/auto/lia/hong/haochuan/NAMD_Git-2021-01-22_Source/cpu_build/Linux-x86_64-g++/namd2
gpath_dir_list=(gpath_[0-9][0-9][0-9][0-9])
for gpath_dir in ${gpath_dir_list[@]};
do
  cd $gpath_dir
  namd_config_list=(committor_gpath_namd_*.namd)
  for namd_config in ${namd_config_list[@]};
  do
    echo "Running $namd_config ..."
    $namd_executable +p1 $namd_config > /dev/null &
  done
  for job in `jobs -p`;
  do
    wait $job
  done
  cd ../
done

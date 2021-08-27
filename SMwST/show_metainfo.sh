#!/bin/sh
string_directory=$1
cd $1/
# get the number of images
num_images=$(grep 'num_images' swarms.conf|awk '{print $3}')
num_iterations=$(grep 'num_iter' swarms.conf|awk '{print $3}')
num_swarm_steps=$(grep 'num_swarm_steps' swarms.conf|awk '{print $3}')
num_equil_steps=$(grep 'num_equil_steps' swarms.conf|awk '{print $3}')
num_total_replicas=$(grep -oP '(?<=replicas )[0-9]+' run.sh)
time_swarm=$(bc <<< "$num_swarm_steps * 0.5")
time_equil=$(bc <<< "$num_equil_steps * 0.5")
num_copies=$(bc <<< "$num_total_replicas / $num_images")
echo "============================================================================"
echo "Simulation $1"
echo "============================================================================"
printf "Number of images: %d\n" $num_images
printf "Number of copies: %d\n" $num_copies
printf "Number of iterations: %d\n" $num_iterations
printf "Simulation time of restrained equilibration: %10.2f fs\n"  $time_swarm
printf "Simulation time of forming swarms of trajectories: %10.2f fs\n"  $time_equil
echo "============================================================================"

#!/bin/sh
png_files=([0-9][0-9].png)
num_files=${#png_files[@]}
num_loops=$(($num_files / 4))
for i in $(seq 1 $num_loops)
do
  file_to_join=""
  for j in $(seq 1 4)
  do
    file_index=$((($i - 1) * 4 + $j))
    file_to_join="$file_to_join $(printf "%02d.png" $file_index)"
  done
  output_name=$(printf "row%02d.png" $i)
  convert +append $file_to_join $output_name
done
# join the rest files
num_rest_files=$(($num_files % 4))
if [[ $num_rest_files -gt 0 ]]
then
  file_to_join=""
  for j in $(seq 1 $num_rest_files)
  do
    file_index=$(($num_loops * 4 + $j))
    file_to_join="$file_to_join $(printf "%02d.png" $file_index)"
  done
  output_name=$(printf "row%02d.png" $(($num_loops + 1)))
  convert +append $file_to_join $output_name
fi

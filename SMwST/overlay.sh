#!/bin/sh
shopt -s nullglob
input_array=(render_[0-9]*.png)
for i in "${!input_array[@]}"
do
    input_image=$(printf "render_%03d.png" $i)
    tmp_image=$(printf "tmp_%03d.png" $i)
    output_image=$(printf "render_overlay_%03d.png" $(($i+1)))
    echo "Processing $input_image..."
    composite -geometry +68 $input_image white_bg.png $tmp_image
    composite template.png $tmp_image $output_image
    rm -f $tmp_image
done

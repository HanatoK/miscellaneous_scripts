#!/bin/sh
shopt -s nullglob
input_array=(render_[0-9][0-9][0-9])
for i in "${!input_array[@]}"
do
    render_filename="${input_array[$i]}"
    output_filename="$render_filename.tga"
    png_filename="$render_filename.png"
    echo "Rendering $render_filename..."
    tachyon -res 1500 1500 $render_filename -o $output_filename
    convert -modulate 100,100,50 -trim $output_filename $png_filename
    rm -f $output_filename
done

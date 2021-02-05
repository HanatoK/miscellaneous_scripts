#!/bin/sh
N_images=35
for i in $(seq 0 $(($N_images/8-1)))
do
  image_row11=$(printf "%02d/dSm_tau_t.png" $((i*8+0)))
  image_row12=$(printf "%02d/dSm_tau_t.png" $((i*8+1)))
  image_row21=$(printf "%02d/dSm_tau_t.png" $((i*8+2)))
  image_row22=$(printf "%02d/dSm_tau_t.png" $((i*8+3)))
  image_row31=$(printf "%02d/dSm_tau_t.png" $((i*8+4)))
  image_row32=$(printf "%02d/dSm_tau_t.png" $((i*8+5)))
  image_row41=$(printf "%02d/dSm_tau_t.png" $((i*8+6)))
  image_row42=$(printf "%02d/dSm_tau_t.png" $((i*8+7)))
  convert +append $image_row11 $image_row12 tmp_row1.png
  convert +append $image_row21 $image_row22 tmp_row2.png
  convert +append $image_row31 $image_row32 tmp_row3.png
  convert +append $image_row41 $image_row42 tmp_row4.png
  convert -append tmp_row1.png tmp_row2.png tmp_row3.png tmp_row4.png page_$i.png
  rm tmp_row[1-4].png
done
remaining_figures=$(($N_images%8))
index=0
for i in $(seq 0 $(($remaining_figures/2-1)))
do
  image_col1=$(printf "%02d/dSm_tau_t.png" $(($N_images/8*8+$i*2+0)))
  image_col2=$(printf "%02d/dSm_tau_t.png" $(($N_images/8*8+$i*2+1)))
  convert +append $image_col1 $image_col2 tmp_row$(($i+1)).png
  index=$(($index+1))
done
echo "Last index: $index"
if (($(($remaining_figures%2)) > 0)); then
  image_last=$(printf "%02d/dSm_tau_t.png" $(($N_images-1)))
  cp $image_last tmp_row$(($index+1)).png
fi
convert -append tmp_row[1-4].png page_$(($N_images/8)).png
rm tmp_row[1-4].png

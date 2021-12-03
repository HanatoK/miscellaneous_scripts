# read colors into a list
set fp [open "colors.dat" "r"]
set file_data [read $fp]
set color_data [split $file_data "\n"]
close $fp

# get a value x from 0 to 1
# return an RGB tuple ranging from 0 to 1
proc get_color_rgb {x} {
  global color_data
  set color_index [expr int(round($x * ([llength $color_data] - 2)))]
  set color_rgba [lindex $color_data $color_index]
  puts "Color is $color_rgba"
  return [lrange $color_rgba 0 2]
}

source plot_path_toc.tcl
set path_files [lsort [glob "path_\[0-9\]\[0-9\]\[0-9\]"]]
# set num_iterations [llength $path_files]
set num_iterations 100
set color_vmd_start 100
graphics top delete all
graphics top materials on
graphics top material Transparent
for {set i 0} {$i < $num_iterations} {incr i} {
  set path_filename [format "path_%03d" $i]
  set color_pos [expr double($i) / ($num_iterations - 1)]
  set color_rgb [get_color_rgb $color_pos]
  color change rgb $color_vmd_start {*}$color_rgb
  plot_path $path_filename $color_vmd_start
  incr color_vmd_start
}
render Tachyon toc.dat

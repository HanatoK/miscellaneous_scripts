source plot_path.tcl
source draw_axis.tcl
set path_files [lsort [glob "path_\[0-9\]\[0-9\]\[0-9\]"]]
set num_iterations [llength $path_files]
for {set i 0} {$i < $num_iterations} {incr i} {
    set path_filename [format "path_%03d" $i]
    set render_filename [format "render_%03d" $i]
    graphics top delete all
    draw_axis
    plot_path $path_filename
    render Tachyon $render_filename
}

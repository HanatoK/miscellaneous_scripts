set image_dirlist [lsort [glob -type d "./image_*"]]
set copy_dirlist [lsort [glob -type d "[lindex $image_dirlist 0]/*"]]
set iteration_dirlist [lsort [glob -type f "[lindex $copy_dirlist 0]/*"]]
set num_images [llength $image_dirlist]
set num_copies [llength $copy_dirlist]
set num_iterations [llength $iteration_dirlist]

for {set i_iteration 0} {$i_iteration < $num_iterations} {incr i_iteration} {
  for {set j_copy 0} {$j_copy < $num_copies} {incr j_copy} {
    set mol_id [mol new ../trialanine.parm7]
    for {set k_image 0} {$k_image < $num_images} {incr k_image} {
      set image_filename [format "image_%02d/copy_%02d/%02d.coor" $k_image $j_copy $i_iteration]
      mol addfile $image_filename $mol_id
    }
    set num_frames [molinfo $mol_id get numframes]
    cv molid $mol_id
    cv configfile cv_dihed.in
    set path_output [open [format "path_copy_%02d_%03d.dat" $j_copy $i_iteration] "w"]
    for {set i 0} {$i < $num_frames} {incr i} {
      cv frame $i
      cv update
      set phi1 [format "%10.4f" [cv colvar phi1 value]]
      set phi2 [format "%10.4f" [cv colvar phi2 value]]
      set phi3 [format "%10.4f" [cv colvar phi3 value]]
      puts $path_output "$phi1 $phi2 $phi3"
    }
    cv delete
    close $path_output
    mol delete $mol_id
  }
}

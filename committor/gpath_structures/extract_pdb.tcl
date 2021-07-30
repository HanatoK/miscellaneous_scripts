# read colvars trajectory
set trajfile "../output_eq/trialanine_gpath_eq_b2+1.colvars.traj"
set fp_trajfile [open $trajfile "r"]
set frame_cv_list [list]
while {[gets $fp_trajfile line] >= 0} {
  set fields [regexp -all -inline {\S+} $line]
  if {[llength $fields] == 0} {
    continue
  }
  if {[lindex $fields 0] == "#"} {
    continue
  }
  if {[lindex $fields 0] == "0"} {
    continue
  }
  lappend frame_cv_list [lindex $fields 1]
}
close $fp_trajfile
puts [format "Read %d CV values from %s." [llength $frame_cv_list] $trajfile]
# read DCD trajectory
set dcdfilename "../output_eq/trialanine_gpath_eq_b2+1.dcd"
set parmfile "../trialanine.parm7"
set mol_id [mol new $parmfile]
mol addfile $dcdfilename waitfor -1 $mol_id
set num_frames [molinfo $mol_id get numframes]
puts [format "Read %d frames from %s." $num_frames $dcdfilename]
set all_atoms [atomselect top all]
set k 0
# filter the frames
set fp_frame_cv_filtered [open "frame_cv.dat" "w"]
for {set i_frame 0} {$i_frame < $num_frames} {incr i_frame 1} {
  set cv_value [lindex $frame_cv_list $i_frame]
  if {[expr abs($cv_value-0.49)] < 1.5e-4} {
    puts $fp_frame_cv_filtered [format "Frame %d, s = %15.10f" $i_frame $cv_value]
    $all_atoms frame $i_frame
    set pdb_filename [format "gpath-%04d.pdb" $k]
    $all_atoms writepdb $pdb_filename
    incr k
  }
}
close $fp_frame_cv_filtered

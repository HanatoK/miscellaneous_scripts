set num_images 25
set namd_psf_file "../../trialanine.parm7"
set mol_id [mol new $namd_psf_file]
for {set i_image 0} {$i_image < $num_images} {incr i_image} {
    set image_filename [format "%02d.coor" $i_image]
    mol addfile $image_filename $mol_id
}
set num_frames [molinfo $mol_id get numframes]
set all_atoms [atomselect top all]
for {set i_frame 0} {$i_frame < $num_frames} {incr i_frame} {
    $all_atoms frame $i_frame
    set pdb_filename [format "string-%02d.pdb" $i_frame]
    $all_atoms writepdb $pdb_filename
}

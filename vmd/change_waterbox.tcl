proc rotation_matrix_from_two_vectors {X Y} {
  set vecX [vecnorm $X]
  set vecY [vecnorm $Y]
  set cos_theta [vecdot $vecX $vecY]
  set vec_norm [vecnorm [veccross $vecX $vecY]]
  set sin_theta [veclength $vec_norm]
  # from wikipedia: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
  set ux [lindex $vec_norm 0]
  set uy [lindex $vec_norm 1]
  set uz [lindex $vec_norm 2]
  set R11 [expr $cos_theta + $ux * $ux * (1.0 - $cos_theta)]
  set R12 [expr $ux * $uy * (1.0 - $cos_theta) - $uz * $sin_theta]
  set R13 [expr $ux * $uz * (1.0 - $cos_theta) + $uy * $sin_theta]
  set R21 [expr $uy * $ux * (1.0 - $cos_theta) + $uz * $sin_theta]
  set R22 [expr $cos_theta + $uy * $uy * (1.0 - $cos_theta)]
  set R23 [expr $uy * $uz * (1.0 - $cos_theta) - $ux * $sin_theta]
  set R31 [expr $uz * $ux * (1.0 - $cos_theta) - $uy * $sin_theta]
  set R32 [expr $uz * $uy * (1.0 - $cos_theta) + $ux * $sin_theta]
  set R33 [expr $cos_theta + $uz * $uz * (1.0 - $cos_theta)]
  set R1 [list $R11 $R12 $R13 0]
  set R2 [list $R21 $R22 $R23 0]
  set R3 [list $R31 $R32 $R33 0]
  set R4 [list 0 0 0 1.0]
  return [list $R1 $R2 $R3 $R4]
}

# load file
set mol1 [mol new "WT.psf"]
mol addfile "thermalization.restart.coor" molid $mol1
# remove water and ions
set nowater [atomselect $mol1 "all not (water or ions)"]
# rotate the molecules:
#   1. find the direction
set ace [atomselect $mol1 "segname ACE"]
set rbd [atomselect $mol1 "segname RBD"]
set com_ace [measure center $ace]
set com_rbd [measure center $rbd]
set direction [vecnorm [vecsub $com_rbd $com_ace]]
#   2. compute the 4x4 matrix
set mat [rotation_matrix_from_two_vectors $direction [list 0.0 0.0 1.0]]
#   3. rotate the selection
$nowater move $mat
# bring to origin
set com [measure center $nowater]
$nowater moveby [vecscale -1.0 $com]
puts "New center: [measure center $nowater]"
# write to file
set psf_filename "nowt.psf"
set pdb_filename "nowt.pdb"
$nowater writepsf $psf_filename
$nowater writepdb $pdb_filename
mol delete $mol1
# solvate the proteins
set solvate_fileprefix "solvated"
solvate $psf_filename $pdb_filename -b 5.0 -minmax {{-60.0 -60.0 -100.0} {60.0 60.0 100.0}} -o $solvate_fileprefix
# ionize
set solvated_psf_filename "${solvate_fileprefix}.psf"
set solvated_pdb_filename "${solvate_fileprefix}.pdb"
autoionize -psf $solvated_psf_filename -pdb $solvated_pdb_filename -sc 0.15 -o "ionized"

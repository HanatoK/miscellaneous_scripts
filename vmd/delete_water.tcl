# select the atoms to be deleted
set wt [atomselect top "water within 8 of protein"]
# build a list of residues from the selected atoms
set remove_residue_list [$wt get residue]
# remove duplicates in the list
set remove_residue_list [lsort -unique $remove_residue_list]
# select all atoms without the above residues
set new_atoms [atomselect top "all not residue ${remove_residue_list}"]
# save PDBF and PSF
$new_atoms writepdb "new_wb.pdb"
$new_atoms writepsf "new_wb.psf"

# Miscellaneous VMD scripts

## GromacsIndexToVMD.tcl

Select atoms in VMD by a GROMACS index file. Usage:

- Open VMD and load your structure.
- Open TK console and run: `source GromacsIndexToVMD.tcl`.
- Select from the GROMACS index file such as `index.ndx`: `SelectionFromIndexFile index.ndx`.
- If your `index.ndx` contains a selection like `[ Protein ]`, then you can access the command in VMD by `$Protein`. For example, you can measure the center of selected atoms by `measure center $Protein`.

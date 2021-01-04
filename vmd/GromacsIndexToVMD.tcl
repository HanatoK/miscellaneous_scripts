proc ParseGromacsIndexFile {filename} {
  set finput [open $filename "r"]
  while { [gets $finput line] >= 0} {
    # trim the leading and ending spaces
    set line [string trim $line]
    if {[string index $line 0] == {[}} {
      # this is a label line with "[ xxx ]"
      set label_name [string trim [string range $line 1 [expr [string length $line] - 2]]]
      dict set result $label_name [list]
    } elseif {[string length $line] > 0} {
      # split and get indices if this line is not empty
      set data [regexp -all -inline {\S+} $line]
      # append indices to the previous list
      set data_prev [dict get $result $label_name]
      lappend data_prev {*}$data
      dict set result $label_name $data_prev
    } else {
      continue
    }
  }
  close $finput
  return $result
}

proc SelectionFromIndexFile {filename} {
  # parse the GROMACS index file
  set d [ParseGromacsIndexFile $filename]
  foreach item [dict keys $d] {
    global $item
    set atom_serial_numbers [dict get $d $item]
    set num_atoms [llength $atom_serial_numbers]
    puts "Selection $item has $num_atoms atoms."
    set $item [atomselect top "serial $atom_serial_numbers"]
    [set $item] global
  }
}

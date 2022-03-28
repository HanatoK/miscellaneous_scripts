# read NAMD restart step from xsc file
proc read_restart_step {xsc_filename} {
  set fp [open ${xsc_filename} "r"]
  set data [split [read $fp] "\n"]
  set step -1
  foreach line $data {
    set s [string trim $line]
    if {[string length $s] > 0} {
      # extract the first character to see if this is a comment
      set first_char [string index $s 0]
      if {[string compare $first_char "#"] == 0} {
        continue
      } else {
        set split_list [regexp -all -inline {\S+} $s]
        set step [lindex $split_list 0]
      }
    }
  }
  return $step
}

proc find_prefix {basename {start_index 0} {old_prefix ""}} {
  # construct the output filename
  set new_prefix $basename
  if {${start_index} > 0} {
    set new_prefix "${basename}+${start_index}"
  }
  set test_filename "${new_prefix}.restart.coor"
  # test file existence
  if {[file exists ${test_filename}] == 1} {
    set new_index [expr ${start_index} + 1]
    # find next index if not exists
    puts "File ${test_filename} exists, test next index ${new_index}"
    return [find_prefix $basename $new_index ${new_prefix}]
  } else {
    return [list ${new_prefix} ${old_prefix}]
  }
}

proc new_prefix {basename {start_index 0}} {
  set p [find_prefix $basename ${start_index} ""]
  return [lindex $p 0]
}

proc old_prefix {basename {start_index 0}} {
  set p [find_prefix $basename ${start_index} ""]
  return [lindex $p 1]
}

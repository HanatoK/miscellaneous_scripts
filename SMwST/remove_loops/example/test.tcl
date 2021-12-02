# load the library
load ../build/libsmwst.so

# helper function for reading files
proc read_pathway_file {path_filename} {
  set result {}
  set fp [open $path_filename r]
  # read file line by line
  while { [gets $fp data] >= 0 } {
    lappend result [regexp -all -inline {\S+} $data]
  }
  close $fp
  return $result
}

# helper function for printing a 2D list
proc print_2d {l} {
  foreach i $l {
    foreach j $i {
      puts -nonewline [format "%15.7f" $j]
    }
    puts -nonewline "\n"
  }
}

# compute drifts of from path1 to path2
proc compute_drifts {path1 path2} {
  set result {}
  foreach image_1 $path1 image_2 $path2 {
    set drift {}
    foreach x1 $image_1 x2 $image_2 {
      lappend drift [expr $x2 - $x1]
    }
    lappend result $drift
  }
  return $result
}

# read image positions from the previous pathway (path_002.txt)
set previous_path [read_pathway_file "path_002.txt"]
set num_images [llength $previous_path]
# initialize N_image zeros for summation
# compute E(X)
set average_drifts [lrepeat $num_images [lrepeat 3 0]]
# compute E(X^2)
set average_drifts2 [lrepeat $num_images [lrepeat 3 0]]
# assume we have run 20 copies
set num_copies 20
for {set i 0} {$i < $num_copies} {incr i} {
  set current_path [read_pathway_file [format "path_copy_%02d_iter_003.txt" $i]]
  # compute the drifts and the squared drifts
  set drifts [compute_drifts $previous_path $current_path]
  puts "============ drifts $i ============"
  print_2d $drifts
  # iterate over the number of images
  for {set j 0} {$j < $num_images} {incr j} {
    # iterate over the number of dimensionality
    for {set k 0} {$k < [llength [lindex $drifts $j]]} {incr k} {
      set current_drift [lindex [lindex $drifts $j] $k]
      set sum [expr [lindex [lindex $average_drifts $j] $k] + $current_drift]
      set sum2 [expr [lindex [lindex $average_drifts2 $j] $k] + $current_drift * $current_drift]
      lset average_drifts "$j $k" $sum
      lset average_drifts2 "$j $k" $sum2
    }
  }
}
# divide average_drifts and average_drifts2 by the number of images
for {set j 0} {$j < $num_images} {incr j} {
  for {set k 0} {$k < [llength [lindex $drifts $j]]} {incr k} {
    lset average_drifts "$j $k" [expr [lindex [lindex $average_drifts $j] $k] / $num_images]
    lset average_drifts2 "$j $k" [expr [lindex [lindex $average_drifts2 $j] $k] / $num_images]
  }
}
# dump the drifts for debugging
puts "========== average drifts =========="
print_2d $average_drifts
puts "========== average drifts2 =========="
print_2d $average_drifts2
# iterate over all images, get randomized new drifts and add them back
puts "========== randomized drifts =========="
# set a factor to decrease the randomness, since as Benoit said, we cannot
# randomize the path all the time, and need to "cool down" the path after a
# few iterations. In a schematic representation, it looks like:
#                         ↑
#        Fully randomized |\
#                         | \
#                         |  \
#                         |   \
#                         |    \
#                         |     \
#                         |      \
#                         |       \
#                         |        \
#                         |         \
#                         |          \
#           Deterministic |           \________
#                         ―――――――――――――――――――――――――>
#                         0    100    200    300
#                                  Iterations
#
# In other words, for a 300-iteration SMwST simulation, I want the first
# 200 iterations use the randomized updating but with a decreasing randomness,
# and the last 100 iterations deterministic. It jusk likes a kind of simulated
# annealing. To implement this, I scale the standard deviation by a factor
# alpha:
#     alpha = max(-3.0 / 2.0 * current_iteration / total_iterations + 1, 0)
#
# assume this is the second iteration
set current_iteration 2
set total_iterations 300
set alpha 0
if {$current_iteration < 200} {
  set alpha [expr -3.0 / 2.0 * $current_iteration / $total_iterations + 1.0]
}
puts "Current alpha factor is  $alpha"
# copy the previous image positions
set new_path [lrange $previous_path 0 end]
for {set j 0} {$j < $num_images} {incr j} {
  for {set k 0} {$k < [llength [lindex $drifts $j]]} {incr k} {
    set mean [lindex [lindex $average_drifts $j] $k]
    # stddev = sqrt(E(X^2) - E(X)^2)
    set stddev [expr sqrt([lindex [lindex $average_drifts2 $j] $k] - $mean * $mean)]
    # scale the standard deviation by alpha
    set scaled_stddev [expr $alpha * $stddev]
    # The command "gaussian" is introduce by libsmwst.so. It can also accept two list.
    # for example, if you run:
    #   set X [gaussian {-1 0 1} {0.5 3.0 2.1}]
    # then you can get a new list X, containing random numbers with mean -1.0, 0.0, 1.0
    # and standard deviation 0.5, 3.0, 2.1, respectively
    set random_drift [gaussian $mean $scaled_stddev]
    puts -nonewline [format "%15.7f" $random_drift]
    # update the new path with drifts
    lset new_path "$j $k" [expr [lindex [lindex $previous_path $j] $k] + $random_drift]
  }
  puts -nonewline "\n"
}
# done. Dump new_path to a file
set outfile [open "path_003.txt" "w"]
puts "=========== Generate new path to path_003.txt ==========="
for {set j 0} {$j < $num_images} {incr j} {
  for {set k 0} {$k < [llength [lindex $drifts $j]]} {incr k} {
    puts -nonewline $outfile [format " %15.7f" [lindex [lindex $new_path $j] $k]]
  }
  puts -nonewline $outfile "\n"
}
close $outfile

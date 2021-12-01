load ../build/libsmwst.so

proc image_distance {i1 i2} {
  set num_elements [llength $i1]
  set sum 0.0
  for {set i 0} {$i < $num_elements} {incr i} {
    set diff [expr [lindex $i1 $i] - [lindex $i2 $i]]
    set sum [expr $sum + $diff * $diff]
  }
  return [expr sqrt($sum)]
}

set new_path [remove_loops_graph {59.7152   -69.7600    60.6074} {34.0876   -71.9973    61.0058}   {  3.1298   -69.0169    59.1743}   {-21.5822   -72.6507    59.1450}   {-51.6666   -72.1258    58.4655}   {-72.8512   -74.2641    56.0041}   {-68.8322   -72.3413    29.2765}   {-69.3674   -74.0722    -0.7440}   {-65.4876   -69.5760   -27.6013}   {-65.4298   -64.3366   -59.3259}   {-61.7858   -62.5816   -75.9464}   {-70.0387   -35.3013   -75.8728}   {-70.1088    -5.0869   -77.6745}   {-71.0188    20.5439   -74.3000}   {-71.3401    49.8733   -75.5665}   {-71.4789    51.9952   -75.7907}   {-71.7494    23.5037   -76.5862}   {-71.3459    -4.6608   -75.2980}   {-71.0885   -33.6515   -79.0025}   {-68.1521   -61.7685   -78.1541}   {-72.6923   -60.2084   -75.9877}   {-71.6801   -30.7226   -77.0749}   {-71.8593    -3.7915   -78.1123}   {-72.6118    24.0721   -78.7021}   {-72.2480    53.0276   -74.9568}   {-55.9298    62.9437   -78.2810}   {-27.1913    59.1967   -80.3898}   {  1.4571    61.2552   -76.6771}   { 28.7506    59.9808   -76.0503}   { 56.8414    59.6899   -77.0058}   { 32.1921    61.0424   -73.8419}   {  3.5253    60.0684   -76.9817}   {-25.5429    60.4825   -79.4977}   {-52.1788    60.2076   -79.1025}   {-70.1944    59.9086   -70.2572}]

set new_path {{60.03051481469165 59.76878325846921 -69.43443060392075} {57.52561728753364 20.421069525692683 -72.81266779503719} {63.15627123721231 -66.45732434109858 -69.18672161366457} {55.58734600107323 -74.08108070359096 -67.2822058006303} {41.29049087707931 -77.43606212191975 -56.32511596695358} {-44.28092969619283 -84.0715120314512 -76.71414723910648} {-78.85063411187653 -87.52659140920188 -69.82194948212863} {-63.91623491981212 -66.61668585176979 -44.512331747247075} {-66.70290378482538 -66.57508528487037 30.014344277202866} {-70.11388300158194 -69.87117652221953 60.218006595617915}}

foreach path_image $new_path {
  foreach elem $path_image {
    puts -nonewline "$elem "
  }
  puts ""
}

# smooth

set new_path [smooth_pathway 0.3 {*}$new_path]
puts "After smoothing"
foreach path_image $new_path {
  foreach elem $path_image {
    puts -nonewline "$elem "
  }
  puts ""
}

set reparam_path [reparametrize [llength $new_path] {*}$new_path]

puts "=====Reparametrized pathway====="
set num_images [llength $reparam_path]
puts "Number of images: $num_images"

set first_image [lindex $reparam_path 0]
puts [format "%15.7f%15.7f%15.7f" [lindex $first_image 0] [lindex $first_image 1] [lindex $first_image 2]]

for {set i 1} {$i < $num_images} {incr i} {
  set previous_image [lindex $reparam_path [expr $i - 1]]
  set current_image [lindex $reparam_path $i]
  set distance [image_distance $previous_image $current_image]
  puts [format "%15.7f%15.7f%15.7f%15.7f" [lindex $current_image 0] [lindex $current_image 1] [lindex $current_image 2] $distance]
}

# moving by reparametrization
if {$num_images == [llength $new_path]} {
  puts "Image drifts after reparametrization:"
  for {set i 0} {$i < $num_images} {incr i} {
    set previous_image [lindex $new_path $i]
    set current_image [lindex $reparam_path $i]
    set distance [image_distance $previous_image $current_image]
    puts [format "%15.7f" $distance]
  }
}

# test Gaussian random numbers
set mean {1.0 2.0 3.0}
set sigma {0.5 0.9 5.0}
set rand_x [gaussian $mean $sigma]
puts $rand_x

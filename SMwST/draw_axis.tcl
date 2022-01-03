proc draw_axis {} {
  graphics top color black
  graphics top line {-180.0 -180.0 -180.0} {-180.0 -180.0 180.0} width 4 style solid
  graphics top line {-180.0 -180.0 -180.0} {-180.0 180.0 -180.0} width 4 style solid
  graphics top line {-180.0 -180.0 -180.0} {180.0 -180.0 -180.0} width 4 style solid
  graphics top line {-180.0 -180.0 180.0} {-180.0 180.0 180.0} width 4 style solid
  graphics top line {-180.0 -180.0 180.0} {180.0 -180.0 180.0} width 4 style solid
  graphics top line {-180.0 180.0 -180.0} {180.0 180.0 -180.0} width 4 style solid
  graphics top line {-180.0 180.0 -180.0} {-180.0 180.0 180.0} width 4 style solid
  graphics top line {180.0 -180.0 -180.0} {180.0 -180.0 180.0} width 4 style solid
  graphics top line {180.0 -180.0 -180.0} {180.0 180.0 -180.0} width 4 style solid
  graphics top line {180.0 180.0 180.0} {-180.0 180.0 180.0} width 4 style solid
  graphics top line {180.0 180.0 180.0} {180.0 -180.0 180.0} width 4 style solid
  graphics top line {180.0 180.0 180.0} {180.0 180.0 -180.0} width 4 style solid
}

proc draw_ticks {} {
  graphics top color black
  # phi_1 ticks
  set firsttime 1
  for {set i -180.0} {$i <= 180.0} {set i [expr $i + 60.0]} {
    graphics top line "$i -180.0 165.0" "$i -180.0 180.0" width 4 style solid
    if {$firsttime < 1} {
      set minor_tick_interval [expr 60.0 / 5.0]
      for {set j [expr $major_tick_previous + $minor_tick_interval]} {$j < $i} {set j [expr $j + $minor_tick_interval]} {
        graphics top line "$j -180.0 172.5" "$j -180.0 180.0" width 2 style solid
      }
    } else {
      set firsttime 0
    }
    set major_tick_previous $i
  }
  # phi_2 ticks
  set firsttime 1
  for {set i -180.0} {$i <= 180.0} {set i [expr $i + 60.0]} {
    graphics top line "165.0 -180.0 $i" "180.0 -180.0 $i" width 4 style solid
    if {$firsttime < 1} {
      set minor_tick_interval [expr 60.0 / 5.0]
      for {set j [expr $major_tick_previous + $minor_tick_interval]} {$j < $i} {set j [expr $j + $minor_tick_interval]} {
        graphics top line "172.5 -180.0 $j" "180.0 -180.0 $j" width 2 style solid
      }
    } else {
      set firsttime 0
    }
    set major_tick_previous $i
  }
  # phi_3 ticks
  set firsttime 1
  for {set i -180.0} {$i <= 180.0} {set i [expr $i + 60.0]} {
    graphics top line "-180.0 $i 180.0" "-165.0 $i 180.0" width 4 style solid
    if {$firsttime < 1} {
      set minor_tick_interval [expr 60.0 / 5.0]
      for {set j [expr $major_tick_previous + $minor_tick_interval]} {$j < $i} {set j [expr $j + $minor_tick_interval]} {
        graphics top line "-180.0 $j 180.0" "-172.5 $j 180.0" width 2 style solid
      }
    } else {
      set firsttime 0
    }
    set major_tick_previous $i
  }
}

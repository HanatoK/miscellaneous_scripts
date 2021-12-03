proc hasPBCEffect {x1 x2 lowerbound upperbound} {
#     puts "x1 = $x1 ; x2 = $x2 \{$lowerbound $upperbound\}"
    set period [expr $upperbound - $lowerbound]
    set half_period [expr 0.5 * $period]
    set dist [expr abs($x1 - $x2)]
#     puts "Dist = $dist ; half_period = $half_period"
    if {$dist <= $half_period} {
        return 0
    } else {
        puts "here is pbc"
        return 1
    }
}

proc drawLinePBC {P1 P2} {
#     puts "P1: $P1"
#     puts "P2: $P2"
    set pbc -1
    set P1m $P1
    set P2m $P2
    set i 0
    foreach x $P1 y $P2 {
#         puts "here"
        if {[hasPBCEffect $x $y -180.0 180.0] > 0} {
#             puts "here"
            if {$x < 0} {
                lset P1m $i -180
            } else {
                lset P1m $i 180
            }
            if {$y < 0} {
                lset P2m $i -180
            } else {
                lset P2m $i 180
            }
            set pbc 1
        }
        incr i
    }
    if {$pbc > -1} {
        puts "PBC line!: $P1 $P1m"
        puts "PBC line!: $P2 $P2m"
        graphics top line $P1 $P1m width 4
        graphics top line $P2m $P2 width 4
    } else {
        graphics top line $P1 $P2 width 4
    }
}

proc plot_path {filename c} {
    set fp [open $filename r]
    set file_data [read $fp]
    set data [split $file_data "\n"]
    set is_first_line 1
    foreach line $data {
        set line [string trim $line]
        if {$line eq ""} {
        } else {
            if {$is_first_line > 0} {
                graphics top color $c
                set is_first_line 0
                set point_x [regexp -inline -all -- {\S+} $line]
                set x1 [lindex $point_x 0]
                set x2 [lindex $point_x 1]
                set x3 [lindex $point_x 2]
            } else {
                graphics top color $c
                set point_y $point_x
                set y1 [lindex $point_y 0]
                set y2 [lindex $point_y 1]
                set y3 [lindex $point_y 2]
                set point_x [regexp -inline -all -- {\S+} $line]
                set x1 [lindex $point_x 0]
                set x2 [lindex $point_x 1]
                set x3 [lindex $point_x 2]
#                 puts "$point_y $point_x"
#                 graphics top line $point_x $point_y width 4
                drawLinePBC $point_x $point_y
            }
#             graphics top color black
            set point_x [regexp -inline -all -- {\S+} $line]
#             graphics top sphere $point_x radius 3
        }
    }
    close $fp
}

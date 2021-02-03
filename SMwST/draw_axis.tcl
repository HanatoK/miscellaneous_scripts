proc draw_axis {} {
    graphics top color black
    graphics top line {-180.0 -180.0 -180.0} {-180.0 -180.0 180.0} width 2 style solid
    graphics top line {-180.0 -180.0 -180.0} {-180.0 180.0 -180.0} width 2 style solid
    graphics top line {-180.0 -180.0 -180.0} {180.0 -180.0 -180.0} width 2 style solid
    graphics top line {-180.0 -180.0 180.0} {-180.0 180.0 180.0} width 2 style solid
    graphics top line {-180.0 -180.0 180.0} {180.0 -180.0 180.0} width 2 style solid
    graphics top line {-180.0 180.0 -180.0} {180.0 180.0 -180.0} width 2 style solid
    graphics top line {-180.0 180.0 -180.0} {-180.0 180.0 180.0} width 2 style solid
    graphics top line {180.0 -180.0 -180.0} {180.0 -180.0 180.0} width 2 style solid
    graphics top line {180.0 -180.0 -180.0} {180.0 180.0 -180.0} width 2 style solid
    graphics top line {180.0 180.0 180.0} {-180.0 180.0 180.0} width 2 style solid
    graphics top line {180.0 180.0 180.0} {180.0 -180.0 180.0} width 2 style solid
    graphics top line {180.0 180.0 180.0} {180.0 180.0 -180.0} width 2 style solid
}

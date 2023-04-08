define(`create_colvars_harmonic_restraint',`harmonics {
    name          restraint_$1
    colvars       $1
    centers       $2
    forceConstant $3
    outputEnergy  on
}

')dnl
create_colvars_harmonic_restraint(colvars, center, constant)
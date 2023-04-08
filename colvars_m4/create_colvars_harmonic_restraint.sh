#!/usr/bin/env sh
# backup old file
mv -f test.colvars test.colvars.bak
m4 -D colvars=dist_1 -D center=0.5 -D constant=10.0 ./create_colvars_harmonic_restraint.m4 >> test.colvars
m4 -D colvars=dist_1 -D center=0.3 -D constant=20.0 ./create_colvars_harmonic_restraint.m4 >> test.colvars

#!/bin/sh
g++ -shared misc.cpp rmsd.cpp -o libmisc.so -DUSE_TCL_STUBS -Wall -fPIC -I/usr/include/eigen3 -I/opt/tcl85/include -L/opt/tcl85/lib64 -ltclstub8.5 -std=c++11 -O2

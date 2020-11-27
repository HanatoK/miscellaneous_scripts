#ifndef MISC_H
#define MISC_H

#include "rmsd.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <tcl.h>

extern "C" {
// minimal_rmsd {{1 2 3} {4 5 6}} {{-2 5 3} {4 1 0}}
int Minimal_rmsd_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                     Tcl_Obj *const objv[]);
int Optimal_rotate_b_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                         Tcl_Obj *const objv[]);
int Optimal_rotation_matrix_btoa_Cmd(ClientData cdata, Tcl_Interp *interp,
                                     int objc, Tcl_Obj *const objv[]);
// compute 3x3 matrix * 3d vector)
void matrix_multiply_vector(double **mat, double *vec, double *res_vec);
// rotate_atoms_by_matrix $matrix $atom_positions
int Rotate_by_matrix_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                         Tcl_Obj *const objv[]);
int Bring_to_center_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                        Tcl_Obj *const objv[]);
int Center_of_geometry_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                           Tcl_Obj *const objv[]);
int Negative_shift_vectors_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                               Tcl_Obj *const objv[]);
int Print_args_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                   Tcl_Obj *const objv[]);
int Vecmult_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                Tcl_Obj *const objv[]);
int Vecsqrt_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                Tcl_Obj *const objv[]);
double Gaussian(double mean, double sigma);
int Gaussian_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                 Tcl_Obj *const objv[]);
int Misc_Init(Tcl_Interp *interp);
}

#endif // MISC_H

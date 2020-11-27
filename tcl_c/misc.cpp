#include "misc.h"

extern "C" {
// minimal_rmsd {{1 2 3} {4 5 6}} {{-2 5 3} {4 1 0}}
int Minimal_rmsd_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                     Tcl_Obj *const objv[]) {
  int atom_set_length[2];
  Tcl_Obj **atom_set_elem[2];
  Tcl_ListObjGetElements(interp, objv[1], &atom_set_length[0],
                         &atom_set_elem[0]);
  Tcl_ListObjGetElements(interp, objv[2], &atom_set_length[1],
                         &atom_set_elem[1]);
  Tcl_Obj ***atom_a_elem =
      (Tcl_Obj ***)malloc(atom_set_length[0] * sizeof(Tcl_Obj **));
  Tcl_Obj ***atom_b_elem =
      (Tcl_Obj ***)malloc(atom_set_length[0] * sizeof(Tcl_Obj **));
  double **atom_set_a =
      (double **)malloc(atom_set_length[0] * sizeof(double *));
  double **atom_set_b =
      (double **)malloc(atom_set_length[0] * sizeof(double *));
  int atom_length; // expected 3
  for (int i = 0; i < atom_set_length[0]; ++i) {
    atom_set_a[i] = (double *)malloc(3 * sizeof(double));
    atom_set_b[i] = (double *)malloc(3 * sizeof(double));
    Tcl_ListObjGetElements(interp, atom_set_elem[0][i], &atom_length,
                           &atom_a_elem[i]);
    Tcl_ListObjGetElements(interp, atom_set_elem[1][i], &atom_length,
                           &atom_b_elem[i]);
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][0], &(atom_set_a[i][0]));
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][1], &(atom_set_a[i][1]));
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][2], &(atom_set_a[i][2]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][0], &(atom_set_b[i][0]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][1], &(atom_set_b[i][1]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][2], &(atom_set_b[i][2]));
  }
  double result = minimalRMSD(atom_set_a, atom_set_b, atom_set_length[0]);
  Tcl_Obj *res = Tcl_NewDoubleObj(result);
  Tcl_SetObjResult(interp, res);
  for (int i = 0; i < atom_set_length[0]; ++i) {
    free(atom_set_a[i]);
    free(atom_set_b[i]);
  }
  free(atom_a_elem);
  free(atom_b_elem);
  free(atom_set_a);
  free(atom_set_b);
  return TCL_OK;
}

int Optimal_rotate_b_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                         Tcl_Obj *const objv[]) {
  int atom_set_length[2];
  Tcl_Obj **atom_set_elem[2];
  Tcl_ListObjGetElements(interp, objv[1], &atom_set_length[0],
                         &atom_set_elem[0]);
  Tcl_ListObjGetElements(interp, objv[2], &atom_set_length[1],
                         &atom_set_elem[1]);
  Tcl_Obj ***atom_a_elem =
      (Tcl_Obj ***)malloc(atom_set_length[0] * sizeof(Tcl_Obj **));
  Tcl_Obj ***atom_b_elem =
      (Tcl_Obj ***)malloc(atom_set_length[0] * sizeof(Tcl_Obj **));
  double **atom_set_a =
      (double **)malloc(atom_set_length[0] * sizeof(double *));
  double **atom_set_b =
      (double **)malloc(atom_set_length[0] * sizeof(double *));
  double **atom_set_b_rotated =
      (double **)malloc(atom_set_length[0] * sizeof(double *));
  int atom_length; // expected 3
  for (int i = 0; i < atom_set_length[0]; ++i) {
    atom_set_a[i] = (double *)malloc(3 * sizeof(double));
    atom_set_b[i] = (double *)malloc(3 * sizeof(double));
    atom_set_b_rotated[i] = (double *)malloc(3 * sizeof(double));
    Tcl_ListObjGetElements(interp, atom_set_elem[0][i], &atom_length,
                           &atom_a_elem[i]);
    Tcl_ListObjGetElements(interp, atom_set_elem[1][i], &atom_length,
                           &atom_b_elem[i]);
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][0], &(atom_set_a[i][0]));
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][1], &(atom_set_a[i][1]));
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][2], &(atom_set_a[i][2]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][0], &(atom_set_b[i][0]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][1], &(atom_set_b[i][1]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][2], &(atom_set_b[i][2]));
  }
  // rotate B
  optimalRotateAtoB(atom_set_b, atom_set_a, atom_set_b_rotated,
                    atom_set_length[0]);
  Tcl_Obj *res = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < atom_set_length[0]; ++i) {
    Tcl_Obj *atom = Tcl_NewListObj(0, NULL);
    Tcl_ListObjAppendElement(interp, atom,
                             Tcl_NewDoubleObj(atom_set_b_rotated[i][0]));
    Tcl_ListObjAppendElement(interp, atom,
                             Tcl_NewDoubleObj(atom_set_b_rotated[i][1]));
    Tcl_ListObjAppendElement(interp, atom,
                             Tcl_NewDoubleObj(atom_set_b_rotated[i][2]));
    Tcl_ListObjAppendElement(interp, res, atom);
  }
  Tcl_SetObjResult(interp, res);
  for (int i = 0; i < atom_set_length[0]; ++i) {
    free(atom_set_a[i]);
    free(atom_set_b[i]);
    free(atom_set_b_rotated[i]);
  }
  free(atom_a_elem);
  free(atom_b_elem);
  free(atom_set_a);
  free(atom_set_b);
  free(atom_set_b_rotated);
  return TCL_OK;
}

int Optimal_rotation_matrix_btoa_Cmd(ClientData cdata, Tcl_Interp *interp,
                                     int objc, Tcl_Obj *const objv[]) {
  int atom_set_length[2];
  Tcl_Obj **atom_set_elem[2];
  Tcl_ListObjGetElements(interp, objv[1], &atom_set_length[0],
                         &atom_set_elem[0]);
  Tcl_ListObjGetElements(interp, objv[2], &atom_set_length[1],
                         &atom_set_elem[1]);
  Tcl_Obj ***atom_a_elem =
      (Tcl_Obj ***)malloc(atom_set_length[0] * sizeof(Tcl_Obj **));
  Tcl_Obj ***atom_b_elem =
      (Tcl_Obj ***)malloc(atom_set_length[0] * sizeof(Tcl_Obj **));
  double **atom_set_a =
      (double **)malloc(atom_set_length[0] * sizeof(double *));
  double **atom_set_b =
      (double **)malloc(atom_set_length[0] * sizeof(double *));
  int atom_length; // expected 3
  for (int i = 0; i < atom_set_length[0]; ++i) {
    atom_set_a[i] = (double *)malloc(3 * sizeof(double));
    atom_set_b[i] = (double *)malloc(3 * sizeof(double));
    Tcl_ListObjGetElements(interp, atom_set_elem[0][i], &atom_length,
                           &atom_a_elem[i]);
    Tcl_ListObjGetElements(interp, atom_set_elem[1][i], &atom_length,
                           &atom_b_elem[i]);
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][0], &(atom_set_a[i][0]));
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][1], &(atom_set_a[i][1]));
    Tcl_GetDoubleFromObj(interp, atom_a_elem[i][2], &(atom_set_a[i][2]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][0], &(atom_set_b[i][0]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][1], &(atom_set_b[i][1]));
    Tcl_GetDoubleFromObj(interp, atom_b_elem[i][2], &(atom_set_b[i][2]));
  }
  double mat[3][3];
  rotationMatrixAtoB(atom_set_b, atom_set_a, &mat, atom_set_length[0]);
  Tcl_Obj *res_mat = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < 3; ++i) {
    Tcl_Obj *atom = Tcl_NewListObj(0, NULL);
    Tcl_ListObjAppendElement(interp, atom, Tcl_NewDoubleObj(mat[i][0]));
    Tcl_ListObjAppendElement(interp, atom, Tcl_NewDoubleObj(mat[i][1]));
    Tcl_ListObjAppendElement(interp, atom, Tcl_NewDoubleObj(mat[i][2]));
    Tcl_ListObjAppendElement(interp, res_mat, atom);
  }
  Tcl_SetObjResult(interp, res_mat);
  for (int i = 0; i < atom_set_length[0]; ++i) {
    free(atom_set_a[i]);
    free(atom_set_b[i]);
  }
  free(atom_a_elem);
  free(atom_b_elem);
  free(atom_set_a);
  free(atom_set_b);
  return TCL_OK;
}

// compute 3x3 matrix * 3d vector)
void matrix_multiply_vector(double **mat, double *vec, double *res_vec) {
  res_vec[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2];
  res_vec[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2];
  res_vec[2] = mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2];
}

// rotate_atoms_by_matrix $matrix $atom_positions
int Rotate_by_matrix_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                         Tcl_Obj *const objv[]) {
  int length_matrix;
  int length_atoms_positions;
  Tcl_Obj **rotatation_matrix_obj;
  Tcl_Obj **atoms_positions_obj;
  Tcl_ListObjGetElements(interp, objv[1], &length_matrix,
                         &rotatation_matrix_obj);
  Tcl_ListObjGetElements(interp, objv[2], &length_atoms_positions,
                         &atoms_positions_obj);
  double **matrix = (double **)malloc(3 * sizeof(double *));
  Tcl_Obj **matrix_elem[3];
  for (int i = 0; i < 3; ++i) {
    int matrix_cols; // expected 3
    matrix[i] = (double *)malloc(3 * sizeof(double));
    Tcl_ListObjGetElements(interp, rotatation_matrix_obj[i], &matrix_cols,
                           &matrix_elem[i]);
    Tcl_GetDoubleFromObj(interp, matrix_elem[i][0], &(matrix[i][0]));
    Tcl_GetDoubleFromObj(interp, matrix_elem[i][1], &(matrix[i][1]));
    Tcl_GetDoubleFromObj(interp, matrix_elem[i][2], &(matrix[i][2]));
  }
  Tcl_Obj ***atom_elem =
      (Tcl_Obj ***)malloc(length_atoms_positions * sizeof(Tcl_Obj **));
  double **atom_set =
      (double **)malloc(length_atoms_positions * sizeof(double *));
  double **atom_set_rotated =
      (double **)malloc(length_atoms_positions * sizeof(double *));
  Tcl_Obj *res = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < length_atoms_positions; ++i) {
    atom_set[i] = (double *)malloc(3 * sizeof(double));
    atom_set_rotated[i] = (double *)malloc(3 * sizeof(double));
    int num_coordinates; // expected 3
    Tcl_ListObjGetElements(interp, atoms_positions_obj[i], &num_coordinates,
                           &atom_elem[i]);
    Tcl_GetDoubleFromObj(interp, atom_elem[i][0], &(atom_set[i][0]));
    Tcl_GetDoubleFromObj(interp, atom_elem[i][1], &(atom_set[i][1]));
    Tcl_GetDoubleFromObj(interp, atom_elem[i][2], &(atom_set[i][2]));
    matrix_multiply_vector(matrix, atom_set[i], atom_set_rotated[i]);
    Tcl_Obj *atom = Tcl_NewListObj(0, NULL);
    Tcl_ListObjAppendElement(interp, atom,
                             Tcl_NewDoubleObj(atom_set_rotated[i][0]));
    Tcl_ListObjAppendElement(interp, atom,
                             Tcl_NewDoubleObj(atom_set_rotated[i][1]));
    Tcl_ListObjAppendElement(interp, atom,
                             Tcl_NewDoubleObj(atom_set_rotated[i][2]));
    Tcl_ListObjAppendElement(interp, res, atom);
  }
  Tcl_SetObjResult(interp, res);
  for (int i = 0; i < 3; ++i) {
    free(matrix[i]);
  }
  free(matrix);
  for (int i = 0; i < length_atoms_positions; ++i) {
    free(atom_set[i]);
    free(atom_set_rotated[i]);
  }
  free(atom_elem);
  free(atom_set);
  free(atom_set_rotated);
  return TCL_OK;
}

int Bring_to_center_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                        Tcl_Obj *const objv[]) {
  if (objc != 2)
    return TCL_ERROR;
  int atom_set_length;
  Tcl_Obj **atom_set_elem;
  Tcl_ListObjGetElements(interp, objv[1], &atom_set_length, &atom_set_elem);
  Tcl_Obj ***atom_elem =
      (Tcl_Obj ***)malloc(atom_set_length * sizeof(Tcl_Obj **));
  double **atom_set = (double **)malloc(atom_set_length * sizeof(double *));
  int atom_elements; // expected 3
  for (int i = 0; i < atom_set_length; ++i) {
    atom_set[i] = (double *)malloc(3 * sizeof(double));
    Tcl_ListObjGetElements(interp, atom_set_elem[i], &atom_elements,
                           &atom_elem[i]);
    if (atom_elements != 3)
      return TCL_ERROR;
    for (int j = 0; j < atom_elements; ++j) {
      Tcl_GetDoubleFromObj(interp, atom_elem[i][j], &(atom_set[i][j]));
    }
  }
  bringToCenter(atom_set, atom_set_length);
  Tcl_Obj *res = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < atom_set_length; ++i) {
    Tcl_Obj *atom = Tcl_NewListObj(0, NULL);
    for (int j = 0; j < atom_elements; ++j) {
      Tcl_ListObjAppendElement(interp, atom, Tcl_NewDoubleObj(atom_set[i][j]));
    }
    Tcl_ListObjAppendElement(interp, res, atom);
  }
  Tcl_SetObjResult(interp, res);
  for (int i = 0; i < atom_set_length; ++i) {
    free(atom_set[i]);
  }
  free(atom_set);
  free(atom_elem);
  return TCL_OK;
}

int Center_of_geometry_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                           Tcl_Obj *const objv[]) {
  if (objc != 2)
    return TCL_ERROR;
  int atom_set_length;
  Tcl_Obj **atom_set_elem;
  Tcl_ListObjGetElements(interp, objv[1], &atom_set_length, &atom_set_elem);
  Tcl_Obj ***atom_elem =
      (Tcl_Obj ***)malloc(atom_set_length * sizeof(Tcl_Obj **));
  double **atom_set = (double **)malloc(atom_set_length * sizeof(double *));
  int atom_elements; // expected 3
  for (int i = 0; i < atom_set_length; ++i) {
    atom_set[i] = (double *)malloc(3 * sizeof(double));
    Tcl_ListObjGetElements(interp, atom_set_elem[i], &atom_elements,
                           &atom_elem[i]);
    if (atom_elements != 3)
      return TCL_ERROR;
    for (int j = 0; j < atom_elements; ++j) {
      Tcl_GetDoubleFromObj(interp, atom_elem[i][j], &(atom_set[i][j]));
    }
  }
  double *center = (double *)malloc(atom_elements * sizeof(double));
  for (int i = 0; i < atom_set_length; ++i) {
    for (int j = 0; j < atom_elements; ++j) {
      center[j] += atom_set[i][j];
    }
  }
  Tcl_Obj *res = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < atom_elements; ++i) {
    center[i] /= double(atom_set_length);
    Tcl_ListObjAppendElement(interp, res, Tcl_NewDoubleObj(center[i]));
  }
  Tcl_SetObjResult(interp, res);
  for (int i = 0; i < atom_set_length; ++i) {
    free(atom_set[i]);
  }
  free(center);
  free(atom_set);
  free(atom_elem);
  return TCL_OK;
}

int Negative_shift_vectors_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                               Tcl_Obj *const objv[]) {
  if (objc != 3)
    return TCL_ERROR;
  // get vector of vectors
  int atom_set_length;
  Tcl_Obj **atom_set_elem;
  Tcl_ListObjGetElements(interp, objv[1], &atom_set_length, &atom_set_elem);
  Tcl_Obj ***atom_elem =
      (Tcl_Obj ***)malloc(atom_set_length * sizeof(Tcl_Obj **));
  double **atom_set = (double **)malloc(atom_set_length * sizeof(double *));
  int atom_elements = 3; // expected 3
  for (int i = 0; i < atom_set_length; ++i) {
    atom_set[i] = (double *)malloc(atom_elements * sizeof(double));
    Tcl_ListObjGetElements(interp, atom_set_elem[i], &atom_elements,
                           &atom_elem[i]);
    if (atom_elements != 3)
      return TCL_ERROR;
    for (int j = 0; j < atom_elements; ++j) {
      Tcl_GetDoubleFromObj(interp, atom_elem[i][j], &(atom_set[i][j]));
    }
  }
  // get the center vector
  int center_length;
  Tcl_Obj **center_obj;
  Tcl_ListObjGetElements(interp, objv[2], &center_length, &center_obj);
  int num_center_elem = 3;
  double center[num_center_elem];
  for (int i = 0; i < center_length; ++i) {
    Tcl_GetDoubleFromObj(interp, center_obj[i], &(center[i]));
  }
  // shift
  Tcl_Obj *res = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < atom_set_length; ++i) {
    Tcl_Obj *atom = Tcl_NewListObj(0, NULL);
    for (int j = 0; j < atom_elements; ++j) {
      atom_set[i][j] -= center[j];
      Tcl_ListObjAppendElement(interp, atom, Tcl_NewDoubleObj(atom_set[i][j]));
    }
    Tcl_ListObjAppendElement(interp, res, atom);
  }
  Tcl_SetObjResult(interp, res);
  // free resources
  for (int i = 0; i < atom_set_length; ++i) {
    free(atom_set[i]);
  }
  free(atom_set);
  free(atom_elem);
  return TCL_OK;
}

int Print_args_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                   Tcl_Obj *const objv[]) {
  Tcl_Channel stdout_channel = Tcl_GetStdChannel(TCL_STDOUT);
  const char *sep = "\n";
  for (int i = 0; i < objc; ++i) {
    Tcl_WriteObj(stdout_channel, objv[i]);
    Tcl_Write(stdout_channel, sep, strlen(sep));
  }
  return TCL_OK;
}

int Vecmult_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                Tcl_Obj *const objv[]) {
  int num_elem;
  Tcl_Obj **data;
  if (objc - 1 < 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"vec1 vec2 ?vec3? ?vec4? ...");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  double *result = (double *)malloc(num_elem * sizeof(double));
  for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
    if (Tcl_GetDoubleFromObj(interp, data[i_elem], &(result[i_elem])) !=
        TCL_OK) {
      free(result);
      return TCL_ERROR;
    }
  }
  int num_elem2;
  for (int i_vec = 2; i_vec < objc; ++i_vec) {
    if (Tcl_ListObjGetElements(interp, objv[i_vec], &num_elem2, &data) !=
        TCL_OK) {
      free(result);
      return TCL_ERROR;
    }
    if (num_elem != num_elem2) {
      Tcl_SetResult(interp,
                    (char *)"vecmult: two vectors don't have the same size",
                    TCL_STATIC);
      free(result);
      return TCL_ERROR;
    }
    for (int i_elem = 0; i_elem < num_elem2; ++i_elem) {
      double df;
      if (Tcl_GetDoubleFromObj(interp, data[i_elem], &df) != TCL_OK) {
        free(result);
        return TCL_ERROR;
      }
      result[i_elem] *= df;
    }
  }
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
    Tcl_ListObjAppendElement(interp, tcl_result,
                             Tcl_NewDoubleObj(result[i_elem]));
  }
  Tcl_SetObjResult(interp, tcl_result);
  free(result);
  return TCL_OK;
}

int Vecsqrt_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                Tcl_Obj *const objv[]) {
  int num_elem;
  Tcl_Obj **data;
  if (objc - 1 != 1) {
    Tcl_WrongNumArgs(interp, 1, objv,
                     (char *)"must provide exactly one argument.");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  double *result = (double *)malloc(num_elem * sizeof(double));
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
    if (Tcl_GetDoubleFromObj(interp, data[i_elem], &(result[i_elem])) !=
        TCL_OK) {
      free(result);
      return TCL_ERROR;
    } else {
      result[i_elem] = sqrt(result[i_elem]);
      Tcl_ListObjAppendElement(interp, tcl_result,
                               Tcl_NewDoubleObj(result[i_elem]));
    }
  }
  Tcl_SetObjResult(interp, tcl_result);
  free(result);
  return TCL_OK;
}

double Gaussian(double mean, double sigma) {
  static int has_saved = 0;
  static double saved = 0.0;
  if (has_saved > 0) {
    double u2 = saved;
    const double x2 = u2 * sigma + mean;
    has_saved = 0;
    return x2;
  } else {
    double u1, u2, w;
    do {
      u1 = (double)(rand()) / (double)(RAND_MAX)*2.0 - 1;
      u2 = (double)(rand()) / (double)(RAND_MAX)*2.0 - 1;
      w = u1 * u1 + u2 * u2;
    } while (w >= 1 || w == 0);
    const double factor = sqrt(-2.0 * log(w) / w);
    u1 *= factor;
    u2 *= factor;
    saved = u2;
    const double x1 = u1 * sigma + mean;
    has_saved = 1;
    return x1;
  }
}

int Gaussian_Cmd(ClientData cdata, Tcl_Interp *interp, int objc,
                 Tcl_Obj *const objv[]) {
  int num_elem_mean;
  int num_elem_sigma;
  Tcl_Obj **tcl_mean;
  Tcl_Obj **tcl_sigma;
  if (objc - 1 != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"must provide two arguments.");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem_mean, &tcl_mean) !=
      TCL_OK) {
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[2], &num_elem_sigma, &tcl_sigma) !=
      TCL_OK) {
    return TCL_ERROR;
  }
  if (num_elem_mean != num_elem_sigma) {
    return TCL_ERROR;
  }
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i_elem = 0; i_elem < num_elem_mean; ++i_elem) {
    double mean, sigma;
    if (Tcl_GetDoubleFromObj(interp, tcl_mean[i_elem], &mean) != TCL_OK) {
      return TCL_ERROR;
    }
    if (Tcl_GetDoubleFromObj(interp, tcl_sigma[i_elem], &sigma) != TCL_OK) {
      return TCL_ERROR;
    }
    Tcl_ListObjAppendElement(interp, tcl_result,
                             Tcl_NewDoubleObj(Gaussian(mean, sigma)));
  }
  Tcl_SetObjResult(interp, tcl_result);
  return TCL_OK;
}

int Misc_Init(Tcl_Interp *interp) {
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
    return TCL_ERROR;
  }
  if (Tcl_PkgProvide(interp, "Misc", "1.0") == TCL_ERROR) {
    return TCL_ERROR;
  }
  Tcl_CreateObjCommand(interp, "minimal_rmsd", Minimal_rmsd_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "optimal_rotate_b", Optimal_rotate_b_Cmd, NULL,
                       NULL);
  Tcl_CreateObjCommand(interp, "optimal_rotation_matrix_btoa",
                       Optimal_rotation_matrix_btoa_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "rotate_atoms_by_matrix", Rotate_by_matrix_Cmd,
                       NULL, NULL);
  Tcl_CreateObjCommand(interp, "bring_to_center", Bring_to_center_Cmd, NULL,
                       NULL);
  Tcl_CreateObjCommand(interp, "center_of_geometry", Center_of_geometry_Cmd,
                       NULL, NULL);
  Tcl_CreateObjCommand(interp, "negative_shift_vectors",
                       Negative_shift_vectors_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "print_args", Print_args_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "vecmult", Vecmult_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "vecsqrt", Vecsqrt_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "gaussian", Gaussian_Cmd, NULL, NULL);
  return TCL_OK;
}
}

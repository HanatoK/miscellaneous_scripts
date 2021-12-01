#include <tcl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <vector>
#include "Pathway.h"
#include "Reparametrization.h"

/*!
 * @brief Removes loops by checking the distance between neighboring images
 *
 * @details
 * This function calculates the distance matrix between any two images on a
 * pathway, and check if the closest image of i is i-1 or i+1. If not, it will
 * remove images between j (the index of the nearest image of image i) and i-1
 * or i+1. It ought to be noted that this is only a C interface to call the C++
 * function remove_loops_simple() to perform the real calculation.
 *
 * @see remove_loops_simple
 *
 * @param[in] images_in 2D array of the input pathway
 * @param[in] num_images_in number of images of the input pathway
 * @param[in] num_element_in dimensionality of the input images
 * @param[out] images_out array of the output pathway
 * @param[out] num_images_out number of images of the output pathway
 * @param[out] num_element_out dimensionality of the output images
 *
 */
void remove_loops_simple_c_interface(
  double** images_in, int num_images_in, int num_element_in,
  double*** images_out, int* num_images_out, int* num_element_out) {
  // copy the C-style 2D array to a C++ vector
  std::vector<Image> pathway(num_images_in);
  for (int i = 0; i < num_images_in; ++i) {
    pathway[i].mImageIndex = i;
    for (int j = 0; j < num_element_in; ++j) {
      pathway[i].mPosition.push_back(images_in[i][j]);
    }
  }
  const std::vector<Image> new_pathway = remove_loops_simple(pathway);
  *num_images_out = (int)(new_pathway.size());
  *num_element_out = (int)(new_pathway[0].mPosition.size());
  *images_out = (double**)malloc((*num_images_out) * sizeof(double*));
  for (int i = 0; i < *num_images_out; ++i) {
    (*images_out)[i] = (double*)malloc((*num_element_out) * sizeof(double));
    for (int j = 0; j < *num_element_out; ++j) {
      (*images_out)[i][j] = new_pathway[i].mPosition[j];
    }
  }
}

/*!
 * @brief Removes loops by traversing a graph with Dijkstra's algorithm
 *
 * @details
 * This function construct a graph from the pathway, using the images as nodes.
 * If the distance between two images is less than a threshold, then these two
 * images belongs to the same node. The threshold is calculated by a factor,
 * namely \p distance_threshold_factor, multiplying with the average distance
 * of the images on the original pathway. After constructing the graph,
 * Dijkstra's algorithm is used to find out the shortest pathway between the
 * original endpoints of the pathway, which should be acyclic. Also, this is a
 * C interface to call the remove_loops_graph C++ function.
 *
 * @see remove_loops_graph
 *
 * @param[in] images_in 2D array of the input pathway
 * @param[in] num_images_in number of images of the input pathway
 * @param[in] num_element_in dimensionality of the input images
 * @param[in] distance_threshold_factor
 *   a factor for determining the overlapping images
 * @param[out] images_out array of the output pathway
 * @param[out] num_images_out number of images of the output pathway
 * @param[out] num_element_out dimensionality of the output images
 *
 */
void remove_loops_graph_c_interface(
  double** images_in, int num_images_in, int num_element_in,
  double distance_threshold_factor,
  double*** images_out, int* num_images_out, int* num_element_out) {
  std::vector<Image> pathway(num_images_in);
  for (int i = 0; i < num_images_in; ++i) {
    pathway[i].mImageIndex = i;
    for (int j = 0; j < num_element_in; ++j) {
      pathway[i].mPosition.push_back(images_in[i][j]);
    }
  }
  const std::vector<Image> new_pathway = remove_loops_graph(pathway, distance_threshold_factor);
  *num_images_out = (int)(new_pathway.size());
  *num_element_out = (int)(new_pathway[0].mPosition.size());
  *images_out = (double**)malloc((*num_images_out) * sizeof(double*));
  for (int i = 0; i < *num_images_out; ++i) {
    (*images_out)[i] = (double*)malloc((*num_element_out) * sizeof(double));
    for (int j = 0; j < *num_element_out; ++j) {
      (*images_out)[i][j] = new_pathway[i].mPosition[j];
    }
  }
}

/*!
 * @brief Interpolate and re-parametrize a pathway
 *
 * @details
 * The pathway is first interpolated by \p num_images_required * 1000 images,
 * and then \p num_images_required nearly equidistant images are picked out to
 * constitute a new pathway (reparametrization). This is a C interface to
 * initialize an instance of the C++ class Reparametrization.
 *
 * @see Reparametrization
 *
 * @param[in] images_in 2D array of the input pathway
 * @param[in] num_images_in number of images of the input pathway
 * @param[in] num_element_in dimensionality of the input images
 * @param[in] num_images_required number of images of the interpolated pathway
 * @param[out] images_out array of the output pathway
 * @param[out] num_images_out number of images of the output pathway
 * @param[out] num_element_out dimensionality of the output images
 */
void reparametrize_c_interface(
  double** images_in, int num_images_in, int num_element_in, int num_images_required,
  double*** images_out, int* num_images_out, int* num_element_out) {
  Matrix mat(num_images_in, num_element_in);
  for (int i = 0; i < num_images_in; ++i) {
    for (int j = 0; j < num_element_in; ++j) {
      mat(i, j) = images_in[i][j];
    }
  }
  const Reparametrization reparam(mat, size_t(num_images_required), 1000);
  const Matrix result = reparam.compute();
  *num_images_out = (int)(result.numRows());
  *num_element_out = (int)(result.numColumns());
  *images_out = (double**)malloc((*num_images_out) * sizeof(double*));
  for (int i = 0; i < *num_images_out; ++i) {
    (*images_out)[i] = (double*)malloc((*num_element_out) * sizeof(double));
    for (int j = 0; j < *num_element_out; ++j) {
      (*images_out)[i][j] = result(i, j);
    }
  }
}

/*!
 * @brief smooth a pathway
 *
 * @details
 * The new position of image i, z'(i), is determined by:
 *   z'(i) = 0.5 * a * z(i-1) + (1-a) * z(i) + 0.5 * a * z(i+1),
 * where the factor a is \p smooth_param.
 *
 * @param[in] images_in 2D array of the input pathway
 * @param[in] num_images_in number of images of the input pathway
 * @param[in] num_element_in dimensionality of the input images
 * @param[in] smooth_param smoothing factor
 * @param[out] images_out array of the output pathway
 * @param[out] num_images_out number of images of the output pathway
 * @param[out] num_element_out dimensionality of the output images
 */
void smooth_c_interface(
  double** images_in, int num_images_in, int num_element_in,
  double smooth_param, double*** images_out, int* num_images_out,
  int* num_element_out) {
  Matrix mat(num_images_in, num_element_in);
  Matrix A(num_images_in, num_images_in);
  for (int i = 0; i < num_images_in; ++i) {
    for (int j = 0; j < num_element_in; ++j) {
      mat(i, j) = images_in[i][j];
    }
    if (i == 0 || i == num_images_in - 1) {
      A(i, i) = 1.0;
    } else {
      A(i, i) = 1.0 - smooth_param;
      A(i, i-1) = 0.5 * smooth_param;
      A(i, i+1) = 0.5 * smooth_param;
    }
  }
  const Matrix result = A * mat;
  *num_images_out = num_images_in;
  *num_element_out = num_element_in;
  *images_out = (double**)malloc((*num_images_out) * sizeof(double*));
  for (int i = 0; i < *num_images_out; ++i) {
    (*images_out)[i] = (double*)malloc((*num_element_out) * sizeof(double));
    for (int j = 0; j < *num_element_out; ++j) {
      (*images_out)[i][j] = result(i, j);
    }
  }
}

extern "C" {

/*!
 * @brief function for testing TCL C interface
 */
int Print_args_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  Tcl_Channel stdout_channel = Tcl_GetStdChannel(TCL_STDOUT);
  const char* sep = "\n";
  for (int i = 0; i < objc; ++i) {
    Tcl_WriteObj(stdout_channel, objv[i]);
    Tcl_Write(stdout_channel, sep, strlen(sep));
  }
  return TCL_OK;
}

/*!
 * @brief implementation of the vecmult command (element-wise multiplication)
 */
int Vecmult_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  int num_elem;
  Tcl_Obj **data;
  if (objc - 1 < 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"vec1 vec2 ?vec3? ?vec4? ...");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  double *result = (double*)malloc(num_elem * sizeof(double));
  for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
    if (Tcl_GetDoubleFromObj(interp, data[i_elem], &(result[i_elem])) != TCL_OK) {
      free(result);
      return TCL_ERROR;
    }
  }
  int num_elem2;
  for (int i_vec = 2; i_vec < objc; ++i_vec) {
    if (Tcl_ListObjGetElements(interp, objv[i_vec], &num_elem2, &data) != TCL_OK) {
      free(result);
      return TCL_ERROR;
    }
    if (num_elem != num_elem2) {
      Tcl_SetResult(interp, (char *) "vecmult: two vectors don't have the same size", TCL_STATIC);
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
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(result[i_elem]));
  }
  Tcl_SetObjResult(interp, tcl_result);
  free(result);
  return TCL_OK;
}

/*!
 * @brief implementation of the vecsqrt command (element-wise square root)
 */
int Vecsqrt_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  int num_elem;
  Tcl_Obj **data;
  if (objc - 1 != 1) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"must provide exactly one argument.");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  double *result = (double*)malloc(num_elem * sizeof(double));
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
    if (Tcl_GetDoubleFromObj(interp, data[i_elem], &(result[i_elem])) != TCL_OK) {
      free(result);
      return TCL_ERROR;
    } else {
      result[i_elem] = sqrt(result[i_elem]);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(result[i_elem]));
    }
  }
  Tcl_SetObjResult(interp, tcl_result);
  free(result);
  return TCL_OK;
}

/*!
 * @brief get a random number that follows the normal distribution
 *
 * @param[in] mean mean of the normal distribution
 * @param[in] sigma standard deviation of the normal distribution
 * @return a random number drawn from the normal distribution with \p mean
 *         as the mean value and \p sigma as the standard deviation
 */
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
      u1 = (double)(rand())/(double)(RAND_MAX) * 2.0 - 1;
      u2 = (double)(rand())/(double)(RAND_MAX) * 2.0 - 1;
      w = u1 * u1 + u2 * u2;
    } while( w >= 1 || w == 0);
    const double factor = sqrt(-2.0 * log(w) / w);
    u1 *= factor;
    u2 *= factor;
    saved = u2;
    const double x1 = u1 * sigma + mean;
    has_saved = 1;
    return x1;
  }
}

/*!
 * @brief TCL-C interface for getting a list of random numbers following
 *        normal distributions
 */
int Gaussian_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  int num_elem_mean;
  int num_elem_sigma;
  Tcl_Obj **tcl_mean;
  Tcl_Obj **tcl_sigma;
  if (objc - 1 != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"must provide two arguments.");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem_mean, &tcl_mean) != TCL_OK) {
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[2], &num_elem_sigma, &tcl_sigma) != TCL_OK) {
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
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(Gaussian(mean, sigma)));
  }
  Tcl_SetObjResult(interp, tcl_result);
  return TCL_OK;
}

/*!
 * @brief TCL-C interface for removing loops by graph
 */
int Remove_loops_graph_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  int num_elem;
  Tcl_Obj **data;
  if (objc - 1 < 4) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"Need more than 3 images");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  const int num_images = objc - 1;
  double **images = (double**)malloc(num_images * sizeof(double*));
  // printf("Number of images: %d\n", num_images);
  int current_elem;
  for (int i = 0; i < num_images; ++i) {
    if (Tcl_ListObjGetElements(interp, objv[i + 1], &current_elem, &data) != TCL_OK) {
      for (int k = 0; k < i; ++k) {
        free(images[k]);
      }
      free(images);
      return TCL_ERROR;
    } else {
      if (current_elem != num_elem) {
        for (int k = 0; k < i; ++k) {
          free(images[k]);
        }
        free(images);
        return TCL_ERROR;
      } else {
        // printf("Number of elements in each image: %d\n", current_elem);
        images[i] = (double*)malloc(num_elem * sizeof(double));
        for (int j = 0; j < num_elem; ++j) {
          if (Tcl_GetDoubleFromObj(interp, data[j], &(images[i][j])) != TCL_OK) {
            for (int k = 0; k < i; ++k) {
              free(images[k]);
            }
            free(images);
            return TCL_ERROR;
          }
        }
      }
    }
  }
  double **images_new = NULL;
  int num_images_new, num_elem_new;
  remove_loops_graph_c_interface(images, num_images, num_elem, 0.7, &images_new, &num_images_new, &num_elem_new);
  // fprintf(stderr, "New pathway: number of images = %d ; dimensions = %d\n", num_images_new, num_elem_new);
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < num_images_new; ++i) {
    Tcl_Obj *tmp_list = Tcl_NewListObj(0, NULL);
    for (int j = 0; j < num_elem_new; ++j) {
      Tcl_ListObjAppendElement(interp, tmp_list, Tcl_NewDoubleObj(images_new[i][j]));
    }
    Tcl_ListObjAppendElement(interp, tcl_result, tmp_list);
  }
  Tcl_SetObjResult(interp, tcl_result);
  /* free the old pathway */
  for (int i = 0; i < num_images; ++i) {
    free(images[i]);
  }
  free(images);
  /* free the new pathway */
  for (int i = 0; i < num_images_new; ++i) {
    free(images_new[i]);
  }
  free(images_new);
  return TCL_OK;
}

/*!
 * @brief TCL-C interface for removing loops by Benoit's approach
 */
int Remove_loops_simple_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  int num_elem;
  Tcl_Obj **data;
  if (objc - 1 < 4) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"Need more than 3 images");
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[1], &num_elem, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  const int num_images = objc - 1;
  double **images = (double**)malloc(num_images * sizeof(double*));
  // printf("Number of images: %d\n", num_images);
  int current_elem;
  for (int i = 0; i < num_images; ++i) {
    if (Tcl_ListObjGetElements(interp, objv[i + 1], &current_elem, &data) != TCL_OK) {
      for (int k = 0; k < i; ++k) {
        free(images[k]);
      }
      free(images);
      return TCL_ERROR;
    } else {
      if (current_elem != num_elem) {
        for (int k = 0; k < i; ++k) {
          free(images[k]);
        }
        free(images);
        return TCL_ERROR;
      } else {
        // printf("Number of elements in each image: %d\n", current_elem);
        images[i] = (double*)malloc(num_elem * sizeof(double));
        for (int j = 0; j < num_elem; ++j) {
          if (Tcl_GetDoubleFromObj(interp, data[j], &(images[i][j])) != TCL_OK) {
            for (int k = 0; k < i; ++k) {
              free(images[k]);
            }
            free(images);
            return TCL_ERROR;
          }
        }
      }
    }
  }
  double **images_new = NULL;
  int num_images_new, num_elem_new;
  remove_loops_simple_c_interface(images, num_images, num_elem, &images_new, &num_images_new, &num_elem_new);
  // fprintf(stderr, "New pathway: number of images = %d ; dimensions = %d\n", num_images_new, num_elem_new);
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < num_images_new; ++i) {
    Tcl_Obj *tmp_list = Tcl_NewListObj(0, NULL);
    for (int j = 0; j < num_elem_new; ++j) {
      Tcl_ListObjAppendElement(interp, tmp_list, Tcl_NewDoubleObj(images_new[i][j]));
    }
    Tcl_ListObjAppendElement(interp, tcl_result, tmp_list);
  }
  Tcl_SetObjResult(interp, tcl_result);
  /* free the old pathway */
  for (int i = 0; i < num_images; ++i) {
    free(images[i]);
  }
  free(images);
  /* free the new pathway */
  for (int i = 0; i < num_images_new; ++i) {
    free(images_new[i]);
  }
  free(images_new);
  return TCL_OK;
}

/*!
 * @brief TCL-C interface for reparametrization
 */
int Reparametrize_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  int num_images_required;
  Tcl_Obj **data;
  if (objc - 1 < 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"num_images vec1 vec2 vec3 ...");
    return TCL_ERROR;
  }
  if (Tcl_GetIntFromObj(interp, objv[1], &num_images_required) != TCL_OK) {
    return TCL_ERROR;
  }
  const int num_images_old = objc - 2;
  double** images_old = (double**)malloc(num_images_old * sizeof(double*));
  int num_elem;
  for (int i = 0; i < num_images_old; ++i) {
    if (Tcl_ListObjGetElements(interp, objv[i+2], &num_elem, &data) != TCL_OK) {
      for (int k = 0; k < i; ++k) {
        free(images_old[k]);
      }
      free(images_old);
      return TCL_ERROR;
    } else {
      images_old[i] = (double*)malloc(num_elem * sizeof(double));
      for (int j = 0; j < num_elem; ++j) {
        if (Tcl_GetDoubleFromObj(interp, data[j], &(images_old[i][j])) != TCL_OK) {
          for (int k = 0; k < i; ++k) {
            free(images_old[k]);
          }
          free(images_old);
          return TCL_ERROR;
        }
      }
    }
  }
  double** images_new = NULL;
  int num_images_new, num_elem_new;
  reparametrize_c_interface(images_old, num_images_old, num_elem, num_images_required, &images_new, &num_images_new, &num_elem_new);
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < num_images_new; ++i) {
    Tcl_Obj *tmp_list = Tcl_NewListObj(0, NULL);
    for (int j = 0; j < num_elem_new; ++j) {
      Tcl_ListObjAppendElement(interp, tmp_list, Tcl_NewDoubleObj(images_new[i][j]));
    }
    Tcl_ListObjAppendElement(interp, tcl_result, tmp_list);
  }
  Tcl_SetObjResult(interp, tcl_result);
  /* free the old pathway */
  for (int i = 0; i < num_images_old; ++i) {
    free(images_old[i]);
  }
  free(images_old);
  /* free the new pathway */
  for (int i = 0; i < num_images_new; ++i) {
    free(images_new[i]);
  }
  free(images_new);
  return TCL_OK;
}

/*!
 * @brief TCL-C interface for smoothing the pathway
 */
int Smooth_pathway_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  double smooth_param;
  Tcl_Obj **data;
  if (objc - 1 < 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"num_images vec1 vec2 vec3 ...");
    return TCL_ERROR;
  }
  if (Tcl_GetDoubleFromObj(interp, objv[1], &smooth_param) != TCL_OK) {
    return TCL_ERROR;
  }
  const int num_images_old = objc - 2;
  double** images_old = (double**)malloc(num_images_old * sizeof(double*));
  int num_elem;
  for (int i = 0; i < num_images_old; ++i) {
    if (Tcl_ListObjGetElements(interp, objv[i+2], &num_elem, &data) != TCL_OK) {
      for (int k = 0; k < i; ++k) {
        free(images_old[k]);
      }
      free(images_old);
      return TCL_ERROR;
    } else {
      images_old[i] = (double*)malloc(num_elem * sizeof(double));
      for (int j = 0; j < num_elem; ++j) {
        if (Tcl_GetDoubleFromObj(interp, data[j], &(images_old[i][j])) != TCL_OK) {
          for (int k = 0; k < i; ++k) {
            free(images_old[k]);
          }
          free(images_old);
          return TCL_ERROR;
        }
      }
    }
  }
  double** images_new = NULL;
  int num_images_new, num_elem_new;
  // smooth_c_interface(images_old, num_images_old, num_elem, num_images_required, &images_new, &num_images_new, &num_elem_new);
  smooth_c_interface(images_old, num_images_old, num_elem, smooth_param, &images_new, &num_images_new, &num_elem_new);
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i = 0; i < num_images_new; ++i) {
    Tcl_Obj *tmp_list = Tcl_NewListObj(0, NULL);
    for (int j = 0; j < num_elem_new; ++j) {
      Tcl_ListObjAppendElement(interp, tmp_list, Tcl_NewDoubleObj(images_new[i][j]));
    }
    Tcl_ListObjAppendElement(interp, tcl_result, tmp_list);
  }
  Tcl_SetObjResult(interp, tcl_result);
  /* free the old pathway */
  for (int i = 0; i < num_images_old; ++i) {
    free(images_old[i]);
  }
  free(images_old);
  /* free the new pathway */
  for (int i = 0; i < num_images_new; ++i) {
    free(images_new[i]);
  }
  free(images_new);
  return TCL_OK;
}

int Smwst_Init(Tcl_Interp *interp) {
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
    return TCL_ERROR;
  }
  if (Tcl_PkgProvide(interp, "Smwst", "1.0") == TCL_ERROR) {
    return TCL_ERROR;
  }
  Tcl_CreateObjCommand(interp, "print_args", Print_args_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "vecmult", Vecmult_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "vecsqrt", Vecsqrt_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "gaussian", Gaussian_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "remove_loops_simple", Remove_loops_simple_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "remove_loops_graph", Remove_loops_graph_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "reparametrize", Reparametrize_Cmd, NULL, NULL);
  Tcl_CreateObjCommand(interp, "smooth_pathway", Smooth_pathway_Cmd, NULL, NULL);
  return TCL_OK;
}

}

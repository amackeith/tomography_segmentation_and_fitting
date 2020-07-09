%module watershedtools_cpp

%{
  #define SWIG_FILE_WITH_INIT  /* To import_array() below */
  #include "watershedtools_cpp.h"
%}

%include "numpy.i"

/* Required for the NumPy C API. Calls to PyArray_SimpleNew
   will result in a segfault if this is omitted */

%init %{
import_array();
%}



%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
{(double *done_arr_interface, int done_arr_x_interface,
    int done_arr_y_interface, int done_arr_z_interface)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
{(double *ret_arr_interface, int ret_arr_x_interface,
    int ret_arr_y_interface, int ret_arr_z_interface)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
{(double *edm_arr_interface, int edm_arr_x_interface,
    int edm_arr_y_interface, int edm_arr_z_interface)}


%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
{(int *output_array, int maxx_,
    int maxy_, int maxz_)}
%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
{(int *in_arr, int maxx,
    int maxy, int maxz)}




%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
 {(double* output_array, int x, int y, int z)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
 {(double* input_array_1, int a, int b, int c)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
 {(double* input_array_2, int d, int e, int f)}


%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
 {(int* in_arr, int xmax, int ymax, int zmax)}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2)
 {(double *moment_of_inertia, int xmoi, int ymoi)}


%apply (int* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
 {(int* in_arr, int xmax, int ymax, int zmax)}

%apply (double* INPLACE_ARRAY1, int DIM1)
 {(double *com, int comlen)}



%include "watershedtools_cpp.h"

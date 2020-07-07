#ifndef ARGOUT_H
#define ARGOUT_H

void simple_add(
        double *output_array, int x, int y, int z,
        double *input_array_1, int a, int b, int c,
        double *input_array_2, int d, int e, int f);

void grow_labels_interface(double *done_arr_interface,
                           int done_arr_x_interface,
                           int done_arr_y_interface,
                           int done_arr_z_interface,
                           double *ret_arr_interface,
                           int ret_arr_x_interface,
                           int ret_arr_y_interface,
                           int ret_arr_z_interface,
                           double *edm_arr_interface,
                           int edm_arr_x_interface,
                           int edm_arr_y_interface,
                           int edm_arr_z_interface);

void hoshen_kopelman3d_interface(int *in_arr,
                                 int maxx, int maxy, int maxz);


void remove_small_labels_interface(int min_vol,
                                   int *in_arr,
                                   int maxx, int maxy, int maxz);


#endif

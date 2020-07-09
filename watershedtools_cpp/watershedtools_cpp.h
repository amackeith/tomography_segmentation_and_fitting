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
                           int edm_arr_z_interface,
                           bool debug);

void hoshen_kopelman3d_interface(int *in_arr,
                                 int maxx, int maxy, int maxz,
                                 bool debug);


void remove_small_labels_interface(int min_vol,
                                   int *in_arr,
                                   int maxx, int maxy, int maxz);

void calculate_moment_of_inertia(int *in_arr,
                           int xmax,
                           int ymax,
                           int zmax,
                           double *moment_of_inertia,
                           int xmoi,
                           int ymoi,
                           double x_center_of_mass,
                           double y_cneter_of_mass,
                           double z_center_of_mass);


void center_of_mass(int *in_arr,
                    int xmax,
                    int ymax,
                    int zmax,
                    double *com,
                    int comlen);

#endif

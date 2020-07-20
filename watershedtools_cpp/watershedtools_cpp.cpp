#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <set>
#include <bits/stdc++.h>
#include <stdexcept>

/*
Author: Arthur MacKeith
Date: July 7, 2020

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"

There is a python verion of this code in the main project directory.
That code is as close to this code as possible. This code is a transcription of
that code into cpp. The python has much better documentation so to understand
the algorithms, (particularly the labeling scheme for grow_labels), I recommend
looking at the explanation of the algorithm in that code.
*/


using namespace std;
class union_find{
public:

    //implementation of union find algorithm
    //sorry it is a bit sparse
    std::vector<int> labels;
    int number_of_labels;

    union_find(){
        labels.push_back(0);
        number_of_labels = 0;
    }

    int uf_find(int x){

        int y = x;
        while (labels[y] != y){
            y = labels[y]; //this makes y the name of a set
        }

        while (labels[x] != x){
            int z = labels[x];
            labels[x] = y; // this collapses the trees
            x = z;
        }

        return y;
    }

    int uf_union(int x, int y){
        //unions two together by pionting theone at the other;
        if (x == 0 || y ==0){
            std::cout << "Major malfunction in hoshen kopelman\n";
            throw std::invalid_argument("Major malfunction in hoshen kopelman");
            return 0;
        }
        int tmp = uf_find(y);
        labels[uf_find(x)] = tmp;

        return tmp;
    }

    int uf_make_set(){
        labels[0] = labels[0] + 1;
        number_of_labels += 1;
        labels.push_back(labels[0]);
        return labels[0];
    }
};


class grow_labels
{
private:
    double *done_arr;
    int done_arr_x;
    int done_arr_y;
    int done_arr_z;
    double *ret_arr;
    int ret_arr_x;
    int ret_arr_y;
    int ret_arr_z;
    double *edm_arr;
    int edm_arr_x;
    int edm_arr_y;
    int edm_arr_z;
    int xmax;
    int ymax;
    int zmax;
    bool debug;
public:

    // constructor just stores arrays.
    grow_labels(double *done_arr_in,
                     int done_arr_x_in, int done_arr_y_in, int done_arr_z_in,
                    double *ret_arr_in,
                    int ret_arr_x_in, int ret_arr_y_in, int ret_arr_z_in,
                    double *edm_arr_in,
                    int edm_arr_x_in, int edm_arr_y_in, int edm_arr_z_in,
                    bool debug_=false):
            done_arr(done_arr_in),
            done_arr_x(done_arr_x_in),
            done_arr_y(done_arr_y_in),
            done_arr_z(done_arr_z_in),
            ret_arr(ret_arr_in),
            ret_arr_x(ret_arr_x_in),
            ret_arr_y(ret_arr_y_in),
            ret_arr_z(ret_arr_z_in),
            edm_arr(edm_arr_in),
            edm_arr_x(edm_arr_x_in),
            edm_arr_y(edm_arr_y_in),
            edm_arr_z(edm_arr_z_in),
            debug(debug_) {
        xmax = done_arr_x;
        ymax = done_arr_y;
        zmax = done_arr_z;
        process();
    }

    // ###################################
    // helper functions
    double ret_arr_at(int x, int y, int z){
        return ret_arr[ymax*zmax*x + zmax*y + z];
    }

    void ret_arr_set(int x, int y, int z, double s){
        ret_arr[ymax*zmax*x + zmax*y + z] = s;
    }

    double edm_arr_at(int x, int y, int z){
        return edm_arr[ymax*zmax*x + zmax*y + z];
    }

    double done_arr_at(int x, int y, int z){
        return done_arr[ymax*zmax*x + zmax*y + z];
    }

    void done_arr_set(int x, int y, int z, double s){
        done_arr[ymax*zmax*x + zmax*y + z] = s;
    }

    bool legit_addr(int x, int y, int z){
        return x >= 0 && y >= 0 && z >= 0 && x < xmax && y < ymax && z < zmax;
    }


    double sum_matrix(double *m){
        long double s = 0;
        for (int x = 0; x < xmax; x++) {
            for (int y = 0; y < ymax; y++) {
                for (int z = 0; z < zmax; z++) {
                    s += m[ymax*zmax*x + zmax*y + z];
                }
            }
        }

        return s;

    }

    // ###################################



    void process(){

        //ret_arr := png.CopyVolume(lbls.Xmax, lbls.Ymax, lbls.Zmax, lbls)
        //ret is copy of lbls
        // done is all zeros at the begining


        //set up done array
        for (int x = 0; x < xmax; x++){
            //cout << x << " of in the set up of " << xmax << "\n";
            for (int y = 0; y < ymax; y++){
                for (int z = 0; z < zmax; z++){
                    if (ret_arr_at(x,y,z) > 0.0 ){
                        done_arr_set(x,y,z,2.0);
                    } else if (edm_arr_at(x,y,z) == 0.0){
                        done_arr_set(x,y,z, 1.0);
                    } else {
                        done_arr_set(x,y,z, 0.0);
                    }
                }
            }
        }

        // FOR DEBUGGING
        /*
        cout << " setting up done array complete\n";
        cout << "sum of done array: " << sum_matrix(done_arr) << "\n";
        cout << "sum of ret array: " << sum_matrix(ret_arr) << "\n";
        cout << "sum of edm array: " << sum_matrix(edm_arr) << "\n";
        */

        //set up the min heap
        //edm height, then the coordinates
        priority_queue<pair<double, vector<int>> > pq;
        //vector<int> coords{ 10, 20, 30 };
        //pq.push(make_pair(10.0, coords));
        //pair<double, vector<int> top = pq.top();

        for (int x = 0; x < xmax; x++){
            if (x % 50 == 0.0 && debug) {
                cout << x << " of  " << xmax << "\n";
            }
            for (int y = 0; y < ymax; y++){
                for (int z = 0; z < zmax; z++){
                    if (done_arr_at(x,y,z) == 2.0) {
                        double key = edm_arr_at(x, y, z);
                        vector<int> coords{x, y, z};
                        pq.push(make_pair(key, coords));
                    }
                }
            }
        }

        //cout << " priority queue set up\n";

        int i,j,k;
        double height;

        while (!pq.empty()){
            if (pq.size() % 50000 == 0 && debug) {
                cout << " pq length is " << pq.size() << "\n";
            }

            pair<double, vector<int>> vox = pq.top();
            pq.pop();

            i = vox.second[0];
            j = vox.second[1];
            k = vox.second[2];
            height = vox.first;

            vector<vector<int>> local_values;
            vector<int> index_adjuster{-1, 0, 1};

            for (auto &di : index_adjuster){
                for (auto &dj : index_adjuster){
                    for (auto &dk : index_adjuster){
                        int sum_abs = abs(di) + abs(dj) + abs(dk);
                        if (sum_abs == 1 &&
                            legit_addr(i+di, j+dj, k+dk)){
                            vector<int> good_addr{i+di, j+dj, k+dk};
                            local_values.push_back(good_addr);
                        } //add all the neighbors to the local values
                    }
                }
            }

            for (int lcl_cnt = 0; lcl_cnt < local_values.size(); lcl_cnt++){
                int x = local_values[lcl_cnt][0];
                int y = local_values[lcl_cnt][1];
                int z = local_values[lcl_cnt][2];

                double edm_height = edm_arr_at(x,y,z);

                if (done_arr_at(x,y,z) == 0.0 && height >= edm_height){
                    ret_arr_set(x,y,z, ret_arr_at(i,j,k));
                    done_arr_set(x,y,z, 2.0);

                    double key = edm_arr_at(x,y,z);
                    vector<int> coords{x,y,z};
                    pq.push(make_pair(key, coords));
                }
            }

            done_arr_set(i, j, k, 1.0);
        }

    }

};



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
                           bool debug){

    grow_labels(done_arr_interface,
            done_arr_x_interface,
            done_arr_y_interface,
            done_arr_z_interface,
            ret_arr_interface,
            ret_arr_x_interface,
            ret_arr_y_interface,
            ret_arr_z_interface,
            edm_arr_interface,
            edm_arr_x_interface,
            edm_arr_y_interface,
            edm_arr_z_interface,
            debug);

}






class hoshen_kopelman3d
{
public:
    int *arr;
    int xmax;
    int ymax;
    int zmax;
    bool debug;

    hoshen_kopelman3d(int *in_arr,
            int maxx, int maxy, int maxz, bool debug_=false): arr(in_arr),
            xmax(maxx), ymax(maxy), zmax(maxz), debug(debug_){
        process();
    }

    bool legit_addr(int x, int y, int z){
        return x >= 0 && y >= 0 && z >= 0 &&
               x < xmax && y < ymax && z < zmax;
    }

    void process() {
        //the expected input is 0 for non-occupied, 1 for occupied
        //invert it so we know anything < 1 has no label on it yet
        for (int x = 0; x < xmax; x++) {
            for (int y = 0; y < ymax; y++) {
                for (int z = 0; z < zmax; z++) {
                    arr[ymax * zmax * x + zmax * y + z] =
                    -1*arr[ymax * zmax * x + zmax * y + z];
                }
            }
        }



        union_find uf = union_find();
        int non_zeros = 0;

        for (int i = 0; i < xmax; i++) {
            if ( i % 50 == 0 && debug){
            cout << " total traversal " << i << " of " << xmax << "\n";
            }
            for (int j = 0; j < ymax; j++) {
                for (int k = 0; k < zmax; k++) {
                    if (arr[ymax * zmax * i + zmax * j + k] == -1) {

                        set<int> nearby_values;
                        //
                        //check "behind" along each axis to see if they have
                        // been marked as part of a segment;
                        vector<int> index_adjuster{1, 0, -1};
                        for (auto &di : index_adjuster) {
                            for (auto &dj : index_adjuster) {
                                for (auto &dk : index_adjuster) {
                                    int sum_abs = abs(di) + abs(dj) + abs(dk);
                                    if (sum_abs == 1 &&
                                        legit_addr(i + di, j + dj, k + dk)) {
                                        // if it is a legit address subject it
                                        // to testing
                                        int neighbor_value =
                                                arr[zmax * ymax * (i + di) +
                                                    zmax * (j + dj) +
                                                    (k + dk)];



                                        if (neighbor_value > 0) {
                                            non_zeros += 1;
                                            /*cout <<
                                            "neighbor value" <<
                                             neighbor_value << "\n";
                                            cout << (i + di) << " " << (j + dj )
                                             << " " <<  (k + dk) << "\n";
                                             */
                                            nearby_values.insert(
                                                                neighbor_value);
                                        }
                                    } //add all the neighbors to nearby values
                                }
                            }
                        }

                        //nearby_values should be a set of the labels touching
                        //this binary 1 voxel.
                        int s = 0; // sum of set labels
                        int max_set = 0;
                        //
                        set<int>::iterator nbv_it = nearby_values.begin();
                        for (nbv_it = nearby_values.begin() ;
                         nbv_it != nearby_values.end(); nbv_it++) {
                            s += *nbv_it;
                            if (*nbv_it > max_set) {
                                max_set = *nbv_it;
                            }
                        }


                        if (!nearby_values.empty()) {

                            for ( nbv_it = nearby_values.begin();
                                  nbv_it != nearby_values.end();
                                  nbv_it++) {
                            }
                        }





                        if (s == 0) { //no labels touching it yet
                            arr[zmax * ymax * i + zmax * j + k] =
                                                            uf.uf_make_set();

                        } else if (nearby_values.size() == 1 && s > 0) {
                            //there is exactly one label touching it
                            arr[zmax * ymax * i + zmax * j + k] = max_set;
                            //max is
                            //arbirary choice, just so this is consistent
                            // accross runs
                        } else {
                            //there are many sets touching it that must be
                            //unioned together
                            for (nbv_it = nearby_values.begin();
                                nbv_it != nearby_values.end();
                                ++nbv_it) {
                                if (max_set != *nbv_it) {
                                    arr[zmax * ymax * i + zmax * j + k] =
                                            uf.uf_union(max_set, *nbv_it);
                                    //make all the sets the same as the max
                                    //this is also arbitrary, but keeps things
                                    //consistent
                                    //not this is different than the python
                                    //implementation but it should get washed
                                    //out if you use better_labels
                                }

                            }
                        }


                    }
                }
            }
        }
        //now we have to fix all the sets that might not be in once piece;

        int s = 0;
        for (int i = 0; i < xmax; i++) {
            if ( i % 100 == 0 && debug){
            cout << " find step " << i << " of " << xmax << "\n";
            }

            for (int j = 0; j < ymax; j++) {
                for (int k = 0; k < zmax; k++) {
                    int val = arr[ymax * zmax * i + zmax * j + k];
                    if (val != 0) {
                        arr[ymax * zmax * i + zmax * j + k]
                        = uf.uf_find(val);
                        s += 1;
                    }

                }
            }
        }


    }


};

void hoshen_kopelman3d_interface(int *in_arr,
                  int maxx, int maxy, int maxz,
                  bool debug){
    hoshen_kopelman3d(in_arr, maxx, maxy, maxz, debug);
}


void simple_add(
        double *output_array, int x, int y, int z,
        double *input_array_1, int a, int b, int c,
        double *input_array_2, int d, int e, int f) {
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                output_array[i * y * z + j * z + k] =
                        input_array_1[i * y * z + j * z + k] +
                        input_array_2[i * y * z + j * z + k];
            }
        }
    }
}


void remove_small_labels_interface(int min_vol,
                                   int *in_arr,
                                   int maxx, int maxy, int maxz){

    //This was not in the original tool kit, it was just something
    //to improve the speed a bit,
    // there is equivlent python code in the section where this is called
    // in watershed_segmentation.py
    map<int, int> label_count;

    for (int i = 0; i < maxx; i++) {
        for (int j = 0; j < maxy; j++) {
            for (int k = 0; k < maxz; k++) {
                int lbl = in_arr[i * maxy * maxz + j * maxz + k];
                label_count[lbl] += 1;
            }
        }
    }

    map<int, int>::iterator lc_it = label_count.begin();
    set<int> labels_to_be_removed;

    for(lc_it = label_count.begin(); lc_it != label_count.end(); ++lc_it){
        if (lc_it->second < min_vol){
            labels_to_be_removed.insert(lc_it->first);
        }
    }

    for (int i = 0; i < maxx; i++) {
        for (int j = 0; j < maxy; j++) {
            for (int k = 0; k < maxz; k++) {
                int lbl = in_arr[i * maxy * maxz + j * maxz + k];
                if (labels_to_be_removed.find(lbl) !=
                labels_to_be_removed.end()){
                    in_arr[i * maxy * maxz + j * maxz + k] = 0;
                }
            }
        }
    }
}



void calculate_moment_of_inertia(int *in_arr,
                           int xmax,
                           int ymax,
                           int zmax,
                           double *moment_of_inertia,
                           int xmoi,
                           int ymoi,
                           double x_center_of_mass,
                           double y_cneter_of_mass,
                           double z_center_of_mass){

    if (xmoi != 3 || ymoi != 3){
        cout << "The moment of inertia matrix is the wrong shape crashing";
        throw std::invalid_argument("Moment of inertia matrix wrong shape");
    }

    ymoi = 3;

    for (int i = 0; i < xmax; i++) {
        for (int j = 0; j < ymax; j++) {
            for (int k = 0; k < zmax; k++) {
                if (in_arr[i * ymax * zmax + j * zmax + k] == 1){
                    double x, y, z;
                    x = i - x_center_of_mass;
                    y = j - y_cneter_of_mass;
                    z = k - z_center_of_mass;

                    moment_of_inertia[3*0 + 0] = moment_of_inertia[3*0 + 0] + pow(y, 2) + pow(z, 2); //I_xx
                    moment_of_inertia[3*1 + 1] = moment_of_inertia[3*1 + 1] + pow(x, 2) + pow(z, 2); // I_yy
                    moment_of_inertia[3*2 + 2] = moment_of_inertia[3*2 + 2] + pow(x, 2) + pow(y, 2); //I_zz
                    moment_of_inertia[3*0 + 1] = moment_of_inertia[3*0 + 1] - x * y; //I_xy
                    moment_of_inertia[3*1 + 2] = moment_of_inertia[3*1 + 2] - y * z; //I_yz
                    moment_of_inertia[3*0 + 2] = moment_of_inertia[3*0 + 2] - x * z; //I_zx

                }
            }
        }
    }
    moment_of_inertia[3*1 + 0] = moment_of_inertia[3*0 + 1]; //I_yx
    moment_of_inertia[3*2 + 1] = moment_of_inertia[3*1 + 2]; //I_zy
    moment_of_inertia[3*2 + 0] = moment_of_inertia[3*0 + 2]; //I_zx
}

void center_of_mass(int *in_arr,
                    int xmax,
                    int ymax,
                    int zmax,
                    double *com,
                    int comlen){
    double x;
    double y;
    double z; //centers of mass in respective coords
    double count;
    count = 0.0;
    x = 0.0;
    y = 0.0;
    z = 0.0;

    for (int i = 0; i < xmax; i++) {
        for (int j = 0; j < ymax; j++) {
            for (int k = 0; k < zmax; k++) {
                if (in_arr[i * ymax * zmax + j * zmax + k] != 0){
                    count += 1;
                    x += i;
                    y += j;
                    z += k;
                }
            }
        }
    }
    x = x / count;
    y = y / count;
    z = z / count;

    com[0] = x;
    com[1] = y;
    com[2] = z;
}



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
using namespace std;
class union_find{
public:

    //implementation of union find algorithm
    //sorry it is a bit sparce
    std::vector<int> labels;
    int number_of_labels;

    union_find(){
        labels.push_back(0);
        number_of_labels = 0;
    }

    int uf_find(int x){
        //cout << "UF FIND START\n";

        int y = x;
        //cout << "SIZE LABELS " << labels.size() << " x " << x <<"\n";
        while (labels[y] != y){
            y = labels[y]; //this makes y the name of a set
        }
        //cout << "UF FIND END 0 \n";

        while (labels[x] != x){
            int z = labels[x];
            labels[x] = y; // this collapses the trees
            x = z;
        }

        //cout << "UF FIND END\n";
        return y;
    }

    int uf_union(int x, int y){
        //unions two together by pionting theone at the other;
        if (x == 0 || y ==0){
            std::cout << "Major malfunction in hoshen kopelman\n";
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
public:
    grow_labels(double *done_arr_in,
                     int done_arr_x_in, int done_arr_y_in, int done_arr_z_in,
                    double *ret_arr_in,
                    int ret_arr_x_in, int ret_arr_y_in, int ret_arr_z_in,
                    double *edm_arr_in,
                    int edm_arr_x_in, int edm_arr_y_in, int edm_arr_z_in):
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
            edm_arr_z(edm_arr_z_in) {
        xmax = done_arr_x;
        ymax = done_arr_y;
        zmax = done_arr_z;
        cout << " process called \n";
        process();
    }


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

        cout << " setting up done array complete\n";
        cout << "sum of done array: " << sum_matrix(done_arr) << "\n";
        cout << "sum of ret array: " << sum_matrix(ret_arr) << "\n";
        cout << "sum of edm array: " << sum_matrix(edm_arr) << "\n";
        //set up the min heap
        //edm height, then the coordinates
        priority_queue<pair<double, vector<int>> > pq;
        //vector<int> coords{ 10, 20, 30 };
        //pq.push(make_pair(10.0, coords));
        //pair<double, vector<int> top = pq.top();

        for (int x = 0; x < xmax; x++){
            if (x % 50 == 0.0) {
                cout << x << " of 2  " << xmax << "\n";
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

        cout << " priority queue set up\n";

        int i,j,k;
        double height;

        while (!pq.empty()){
            if (pq.size() % 10000 == 0) {
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
                           int edm_arr_z_interface){

    cout << " grow labels interface called\n";
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
            edm_arr_z_interface);

}



class better_labels {
    // the python version of this with numpy is actually faster
    // but I only found out after I wrote this, leaving this here for
    // completeness
    int *arr;
    int *output_arr;
    int xmax;
    int ymax;
    int zmax;
public:
    better_labels(int *processed_array,
                int maxx, int maxy, int maxz,
                int *output_array):
                arr(processed_array), output_arr(output_array),
                xmax(maxx), ymax(maxy), zmax(maxz){
        process();
    }

    void process(){
        int next_label = 0;

        set<int> lbls;
        set<int>::iterator lbls_itter;

        for (int x = 0; x < xmax; x++) {
            for (int y = 0; y < ymax; y++) {
                for (int z = 0; z < zmax; z++) {
                    lbls.insert(arr[ymax*zmax*x + zmax*y + z]);
                }
            }
        }

        for (lbls_itter = lbls.begin(); lbls_itter != lbls.end(); ++lbls_itter){
            if (*lbls_itter != 0){

                next_label += 1;

                for (int x = 0; x < xmax; x++) {
                    for (int y = 0; y < ymax; y++) {
                        for (int z = 0; z < zmax; z++) {
                            if (arr[ymax*zmax*x + zmax*y + z] == *lbls_itter){
                                output_arr[ymax*zmax*x + zmax*y + z] =
                                        next_label;
                            }
                        }
                    }
                }

            }


        }
    }
};


void better_labels_interface(int *in_arr,
                           int maxx, int maxy, int maxz,
                           int *output_array,
                           int maxx_, int maxy_, int maxz_){
    better_labels(in_arr,
                  maxx, maxy, maxz,
                  output_array);
}



class hoshen_kopelman3d
{
public:
    int *arr;
    int xmax;
    int ymax;
    int zmax;

    hoshen_kopelman3d(int *in_arr,
            int maxx, int maxy, int maxz): arr(in_arr),
            xmax(maxx), ymax(maxy), zmax(maxz){
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
            if ( i % 50 == 0){
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
                                            //cout << "neighbor value" << neighbor_value << "\n";
                                            //cout << (i + di) << " " << (j + dj ) << " " <<  (k + dk) << "\n";

                                            nearby_values.insert(neighbor_value);
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
                        for (nbv_it = nearby_values.begin() ; nbv_it != nearby_values.end(); nbv_it++) {
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
                            arr[zmax * ymax * i + zmax * j + k] = uf.uf_make_set();

                        } else if (nearby_values.size() == 1 && s > 0) {
                            //there is exactly one label touching it
                            arr[zmax * ymax * i + zmax * j + k] = max_set; //max is
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
        cout << "SUM OF THINGS " << non_zeros << "\n";
        //now we have to fix all the sets that might not be in once piece;

        int s = 0;
        cout << "it is the find step\n";
        for (int i = 0; i < xmax; i++) {
            if ( i % 50 == 0){
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




        cout << "PRINT" << uf.number_of_labels <<  "\n";
        long int sum = 0 ;
        for (int i = 0; i < xmax; i++) {
            for (int j = 0; j < ymax; j++) {
                for(int k = 0; k < zmax; k++){
                    //
                    if (arr[ymax * zmax * i + zmax * j + k] > 10 || -10 > arr[ymax * zmax * i + zmax * j + k]){
                    // cout << arr[ymax * zmax * i + zmax * j + k] << " " << i << " " << j << " " << k << "\n";
                    }
                    sum += arr[ymax * zmax * i + zmax * j + k];
                }

                }

        }
        cout << "PRINT2\n";


        cout << "sum is " << sum << "\n";

    }


};

void hoshen_kopelman3d_interface(int *in_arr,
                  int maxx, int maxy, int maxz){
    hoshen_kopelman3d(in_arr, maxx, maxy, maxz);
    cout << "interface returns\n";
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
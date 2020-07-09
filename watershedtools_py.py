'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"


There is a cpp version of this code in the watershedtools_cpp folder of this project
That code runs an order of magnitude faster so if you can get it to work it is
well worth it.

This python code is as close as I could make it, so hopefully if the cpp is
difficult to understand the python code can be a guide to what is going on.
=========
'''

import numpy as np
import heapq


# makes the labels count from 1 to N rather than having them all over the map
def better_labels(processed_array):
    # relabels an array (skipping 0) so that labels begin at 1
    # and go up to the number of labels
    arr = np.int32(processed_array.copy())
    
    lbls = arr.copy()
    lbls = lbls.flatten()
    lbls = list(set(lbls))
    lbls.sort()
    
    out_arr = np.zeros(arr.shape)
    
    new_label = 0
    
    for lbl in lbls:
        if lbl == 0:
            continue
        
        new_label = new_label + 1
        out_arr[arr == lbl] = new_label
    
    return np.int32(out_arr)


# these are the tools of union find,
# I implemented this following
# https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
# as a guide
labels = [0, ]
number_of_labels = 0


def uf_find(x):
    x = int(x)
    y = x
    while labels[y] != y:
        y = labels[y]  # this makes y the name of a set
    
    while labels[x] != x:
        z = labels[x]
        labels[x] = y  # this collapses trees
        x = z
    return y


def uf_union(x, y):
    # unions two together by pointing the one at the other
    if x == 0 or y == 0:
        print("major wrong in hoshen")
        exit()
    
    x = int(x);
    y = int(y)
    tmp = uf_find(y)
    labels[uf_find(x)] = tmp
    # print("union", x,y)
    return tmp


def uf_make_set():
    # adds a new label to the label array
    labels[0] = labels[0] + 1
    global number_of_labels
    number_of_labels = number_of_labels + 1
    labels.append(labels[0])
    return labels[0]


def hoshen_kopelmann_verify_3d(array):
    shp = array.shape
    for i in range(shp[0]):
        
        for j in range(shp[1]):
            for k in range(shp[2]):
                if array[i, j, k] > 0:
                    
                    nearby_values = []  # this will contain the nearby sets
                    
                    # north and west and back
                    # these are actually the only ones we care about bc
                    # south and east and forward have not been assigned yet
                    # this is because of how we are walking out into the array
                    # as the leading corner
                    if i != 0 and j != 0 and k != 0:
                        # west
                        nearby_values.append(array[i - 1, j, k])
                        # north
                        nearby_values.append(array[i, j - 1, k])
                        # back
                        nearby_values.append(array[i, j, k - 1])
                    elif i != 0 and j == 0 and k != 0:
                        # west
                        nearby_values.append(array[i - 1, j, k])
                        # no north possible
                        # back
                        nearby_values.append(array[i, j, k - 1])
                    
                    elif i == 0 and j != 0 and k != 0:
                        # no west possible
                        # north
                        nearby_values.append(array[i, j - 1, k])
                        # back
                        nearby_values.append(array[i, j, k - 1])
                    elif i != 0 and j != 0 and k == 0:
                        # west
                        nearby_values.append(array[i - 1, j, k])
                        # north
                        nearby_values.append(array[i, j - 1, k])
                        # no back possible
                    elif i != 0 and j == 0 and k == 0:
                        # only west
                        nearby_values.append(array[i - 1, j, k])
                    elif i == 0 and j != 0 and k == 0:
                        # only north
                        nearby_values.append(array[i, j - 1, k])
                    elif i == 0 and j == 0 and k != 0:
                        # only back
                        nearby_values.append(array[i, j, k - 1])
                    
                    else:
                        t = (i == 0 and j == 0 and k == 0)
                        if not t:
                            print(nearby_values)
                            print("major malfunction1")
                    
                    # east, south, forward (this is so the second pass works)
                    i_max = shp[0] - 1
                    j_max = shp[1] - 1
                    k_max = shp[2] - 1
                    
                    if i != i_max and j != j_max and k != k_max:
                        # east
                        nearby_values.append(array[i + 1, j, k])
                        # south
                        nearby_values.append(array[i, j + 1, k])
                        # forward
                        nearby_values.append(array[i, j, k + 1])
                    elif i != i_max and j == j_max and k != k_max:
                        # east
                        nearby_values.append(array[i + 1, j, k])
                        # no south possible
                        # forward
                        nearby_values.append(array[i, j, k + 1])
                    elif i == i_max and j != j_max and k != k_max:
                        # no east possible
                        # south
                        nearby_values.append(array[i, j + 1, k])
                        # forward
                        nearby_values.append(array[i, j, k + 1])
                    elif i != i_max and j != j_max and k == k_max:
                        # east
                        nearby_values.append(array[i + 1, j, k])
                        # south
                        nearby_values.append(array[i, j + 1, k])
                        # no forward possible
                    elif i != i_max and j == j_max and k == k_max:
                        # only east
                        nearby_values.append(array[i + 1, j, k])
                    elif i == i_max and j != j_max and k == k_max:
                        # only south
                        nearby_values.append(array[i, j + 1, k])
                    elif i == i_max and j == j_max and k != k_max:
                        # only forward
                        nearby_values.append(array[i, j, k + 1])
                    else:
                        t = (i == i_max and j == j_max and k == k_max)
                        if not t:
                            print(nearby_values, i, j, k, t)
                            print("major malfunction2")
                    
                    nearby_values = filter(lambda item: item > 0, nearby_values)
                    nearby_values = list(set(nearby_values))
                    for item in nearby_values:
                        if array[i, j, k] == item:
                            pass
                        else:
                            print("hoshen_kopelmann3d_union_find failed",
                                  i, j, k)


def hoshen_kopelmann3d(original_input_int_array):
    # the expected input is 0 for non-occupied, 1 for occupied
    input_int_array = np.int32(original_input_int_array.copy())
    shp = input_int_array.shape
    input_int_array = -1 * input_int_array
    # change it to -1 and 0 so we know anything < 1 has no label on it yets
    for i in range(shp[0]):
        print("total traversal", i, "of", shp[0])
        for j in range(shp[1]):
            for k in range(shp[2]):
                if input_int_array[i, j, k] == -1:
                    # this will contain the nearby sets names
                    nearby_values = []
                    
                    # north and west and back
                    # these are actually the only ones we care about bc
                    # south and east and forward have not been assigned yet
                    # this is because of how we are walking out into the array
                    # as the leading corner
                    if i != 0 and j != 0 and k != 0:
                        # west
                        nearby_values.append(input_int_array[i - 1, j, k])
                        # north
                        nearby_values.append(input_int_array[i, j - 1, k])
                        # back
                        nearby_values.append(input_int_array[i, j, k - 1])
                    elif i != 0 and j == 0 and k != 0:
                        # west
                        nearby_values.append(input_int_array[i - 1, j, k])
                        # no north possible
                        # back
                        nearby_values.append(input_int_array[i, j, k - 1])
                    
                    elif i == 0 and j != 0 and k != 0:
                        # no west possible
                        # north
                        nearby_values.append(input_int_array[i, j - 1, k])
                        # back
                        nearby_values.append(input_int_array[i, j, k - 1])
                    elif i != 0 and j != 0 and k == 0:
                        # west
                        nearby_values.append(input_int_array[i - 1, j, k])
                        # north
                        nearby_values.append(input_int_array[i, j - 1, k])
                        # no back possible
                    elif i != 0 and j == 0 and k == 0:
                        # only west
                        nearby_values.append(input_int_array[i - 1, j, k])
                    elif i == 0 and j != 0 and k == 0:
                        # only north
                        nearby_values.append(input_int_array[i, j - 1, k])
                    elif i == 0 and j == 0 and k != 0:
                        # only back
                        nearby_values.append(input_int_array[i, j, k - 1])
                    
                    else:
                        t = (i == 0 and j == 0 and k == 0)
                        if not t:
                            print(nearby_values)
                            print("major malfunction1")
                    
                    # remove empty space and duplicate values
                    nearby_values = filter(lambda item: item > 0, nearby_values)
                    nearby_values = list(set(nearby_values))
                    nearby_values.sort()
                    nearby_values.reverse()
                    s = sum(nearby_values)
                    if s == 0:  # no labels touching it yet
                        input_int_array[i, j, k] = uf_make_set()
                    elif (len(nearby_values) == 1 and s > 0):
                        # there is only one non-zero label touchign
                        input_int_array[i, j, k] = nearby_values[0]
                        #only nonzero value
                    else:
                        # there are multiple values and they must be unioned
                        # (this is arbitrary) but consistent with the cpp code
                        for lbl in nearby_values[1:]:
                            #nearby_values[0] is max nearby values
                            input_int_array[i, j, k] = \
                                uf_union(nearby_values[0], lbl)
    
    for x in range(shp[0]):
        print("find step", x, "of", shp[0])
        # print("2",x)
        for y in range(shp[1]):
            for z in range(shp[2]):
                if input_int_array[x, y, z] != 0:
                    input_int_array[x, y, z] = uf_find(input_int_array[x, y, z])
    
    #return the globals to original states
    global labels, number_of_labels
    labels = [0, ]
    number_of_labels = 0
    # input_int_array = better_labels(input_int_array)
    return input_int_array




# this function grows the labels given an edm array and a seed array
# it expects two same shape arrays

# ### edm -> Euclidian distance map.
# in analogy to real watersheds, lbl_arr_in are the drains, and edm_arr_in
# represents the topology of the land (in this case though the "water" would
# flow up the gradient (and the drains are at the local maxima of the edm
# array.
def grow_labels(edm_arr_in, lbl_arr_in):
    lbl_arr = lbl_arr_in.copy()  # array with the labels that will be grown
    edm_arr = edm_arr_in.copy()  # array with "height map" of volume
    if not (np.array_equal(edm_arr.shape, lbl_arr.shape)
            and len(edm_arr.shape) == 3):
        print("grow_labels in hoshen_kopelmann.py passed two "
              "different size arrays or not length 3")
    
    # this is the array that we will update as we go,
    # so that the lbl_arr is not changed in place
    
    ret_arr = lbl_arr.copy()
    
    # this will be the guide that tells us what we have done so far
    # and what is still left to do.
    # dictionary: done =1, unvisited but labeled: 2, unvisited and unlabeled:0
    # done contains both all the background as well as the visited and labeled.
    # this algorithm will have some unvisited and unlabled on voxels left at the
    # end, these are the ones
    # that were part of lbls that are small and unconnected
    # these will be the ones that are unconnected to anything and will be
    # thrown out. ie set to done if
    # there are no more labled but unvisited
    
    # set everything that is zero in edm to 1 (done) in done array
    done_arr = np.zeros(edm_arr.shape)
    done_arr[edm_arr == 0] = 1
    # set all labeled voxels to unvisited, labled (2),
    # the rest (unvisited unlabled) being 0
    done_arr[ret_arr > 0] = 2
    
    # now make a heap of unvisited and labled voxels
    unvisited_labled = []
    
    # we will visit voxels in decending order of height to we grow
    # the lbls out from the center of each lentil
    
    shp = done_arr.shape
    for i in np.arange(done_arr.shape[0]):
        for j in np.arange(done_arr.shape[1]):
            for k in np.arange(done_arr.shape[2]):
                if done_arr[i, j, k] == 2:
                    # ie if this location is unvisited but labeled
                    height = edm_arr[i, j, k]
                    loc = [i, j, k]
                    # make the height negative bc this is a min heap
                    item = (-height, loc)
                    heapq.heappush(unvisited_labled, item)
    
    # sort list and set it so it gets starts with the largest heights
    print(su)
    
    while len(unvisited_labled) > 0:
        if len(unvisited_labled) % 10000 == 0:
            print("grow labels", len(unvisited_labled), "left to go")
        unv_lbl = heapq.heappop(unvisited_labled)
        # now we will properly update neighbors
        # look at the 6 close by voxels, if any are
        # unvisited+unlabled and have a lower height
        # assign them to the label of unv_lbl, otherwise leave them alone
        i, j, k = unv_lbl[1]
        height = -1 * unv_lbl[0]  # remember to make it positive again
        
        nearby_values = []
        # north and west and back
        # sorry about this it is very verbose way to write a filter.
        # in short it checks the three voxesl behind the leading corner
        # of the way this is itterating. and then in front
        if i != 0 and j != 0 and k != 0:
            # west
            nearby_values.append([i - 1, j, k])
            # north
            nearby_values.append([i, j - 1, k])
            # back
            nearby_values.append([i, j, k - 1])
        elif i != 0 and j == 0 and k != 0:
            # west
            nearby_values.append([i - 1, j, k])
            # no north possible
            # back
            nearby_values.append([i, j, k - 1])
        
        elif i == 0 and j != 0 and k != 0:
            # no west possible
            # north
            nearby_values.append([i, j - 1, k])
            # back
            nearby_values.append([i, j, k - 1])
        elif i != 0 and j != 0 and k == 0:
            # west
            nearby_values.append([i - 1, j, k])
            # north
            nearby_values.append([i, j - 1, k])
            # no back possible
        elif i != 0 and j == 0 and k == 0:
            # only west
            nearby_values.append([i - 1, j, k])
        elif i == 0 and j != 0 and k == 0:
            # only north
            nearby_values.append([i, j - 1, k])
        elif i == 0 and j == 0 and k != 0:
            # only back
            nearby_values.append([i, j, k - 1])
        
        else:
            t = (i == 0 and j == 0 and k == 0)
            if not t:
                print(nearby_values)
                print("major malfunction1")
        
        # east, south, forward (this is so the second pass works)
        i_max = shp[0] - 1
        j_max = shp[1] - 1
        k_max = shp[2] - 1
        # print(i_max,j_max,k_max)
        
        if i != i_max and j != j_max and k != k_max:
            # east
            nearby_values.append([i + 1, j, k])
            # south
            nearby_values.append([i, j + 1, k])
            # forward
            nearby_values.append([i, j, k + 1])
        elif i != i_max and j == j_max and k != k_max:
            # east
            nearby_values.append([i + 1, j, k])
            # no south possible
            # forward
            nearby_values.append([i, j, k + 1])
        elif i == i_max and j != j_max and k != k_max:
            # no east possible
            # south
            nearby_values.append([i, j + 1, k])
            # forward
            nearby_values.append([i, j, k + 1])
        elif i != i_max and j != j_max and k == k_max:
            # east
            nearby_values.append([i + 1, j, k])
            # south
            nearby_values.append([i, j + 1, k])
            # no forward possible
        elif i != i_max and j == j_max and k == k_max:
            # only east
            nearby_values.append([i + 1, j, k])
        elif i == i_max and j != j_max and k == k_max:
            # only south
            nearby_values.append([i, j + 1, k])
        elif i == i_max and j == j_max and k != k_max:
            # only forward
            nearby_values.append([i, j, k + 1])
        else:
            t = (i == i_max and j == j_max and k == k_max)
            if not t:
                print(nearby_values, i, j, k, t)
                print("major malfunction2")
        
        for nearby in nearby_values:
            x, y, z = nearby
            
            if done_arr[x, y, z] == 0 and height >= edm_arr[x, y, z]:
                # if it is unlabeled and unvisited, and it is not uphill
                ret_arr[x, y, z] = ret_arr[i, j, k]
                done_arr[x, y, z] = 2
                
                htxyz = edm_arr[x, y, z]
                locxyz = [x, y, z]
                item = (-htxyz, locxyz)
                heapq.heappush(unvisited_labled, item)
        
        # say I am done looking at i,j,k
        done_arr[i, j, k] = 1
        # print(len(unvisited_labled), height)
    
    ret_arr = better_labels(ret_arr)
    return ret_arr




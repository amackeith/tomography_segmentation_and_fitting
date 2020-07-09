'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"
'''

import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import sys
from moment_of_inertia_orientation import find_principal_axis, get_phi_and_theta_from_principle_axis
import skimage.feature
from multiprocessing import Process, Queue
import multiprocessing as mp
import time
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import sys


cpp_enabled = True
try:
    import watershedtools_cpp.watershedtools_cpp as watershedtools_cpp
except:
    print("CPP code not enabled")
    print("to build and enable it run ./build in the watershedtools_cpp folder")
    print("Proceeding using SLOWER python code")
    cpp_enabled = False
    
    
def get_loc_and_angle(full_fitt, lbl_lst, que, core_num,
                      padding, debug, oblate):
    com_angle_list = []
    num = len(lbl_lst)
    cnt = 0
    
    for particle_index in lbl_lst:
        mask = np.int32((full_fitt == particle_index) * 1)
        mask_vol = np.sum(mask)
        np_where_info = np.where(mask == 1)
        x_max = max(np_where_info[0])
        x_min = min(np_where_info[0])
        y_max = max(np_where_info[1])
        y_min = min(np_where_info[1])
        z_max = max(np_where_info[2])
        z_min = min(np_where_info[2])
        
        x_avg = (x_min + x_max) // 2
        y_avg = (y_min + y_max) // 2
        z_avg = (z_min + z_max) // 2
        
        arr = mask[x_avg - padding:x_avg + padding,
                   y_avg - padding:y_avg + padding,
                   z_avg - padding:z_avg + padding]
        arr = np.ascontiguousarray(arr, dtype=np.int32)
        
        val, vec, com = find_principal_axis(arr)
        
        if oblate: #particles are oblate
            index_of_max_eigen_val = np.argmax(val)
        else: #prolate particle
            index_of_max_eigen_val = np.argmin(val)
            
        principle_axis = vec[:, index_of_max_eigen_val]
        phi, theta = get_phi_and_theta_from_principle_axis(principle_axis)
        
        correctd_com = com - padding  # to correct for nesting it in the empty array
        
        lentil = [correctd_com, np.array([phi, theta, 0]), principle_axis, mask_vol, particle_index, (val, vec)]
        com_angle_list.append(lentil)
        
        if debug:
            print("core_num", core_num, "number", cnt, "of", num, "mask vol", mask_vol)
        cnt = cnt + 1
    
    que.put(com_angle_list)



class orientations_from_moment_of_inertia:
    
    def __init__(self, fname, segmented_volume, padding=30,
                oblate=True, debug=False, force_sequential=False,
                 specify_max_cores=0):
        self.fname = fname[:-4]
        self.volume = np.int32(segmented_volume)
        self.padding = padding
        self.oblate = oblate
        self.debug = debug
        self.force_sequential = force_sequential
        self.specify_max_cores = specify_max_cores
        
        # inset the volume with #padding voxels in each direction
        shp = np.array(self.volume.shape) + 2 * padding
        nest = np.zeros(shp)
        nest[padding:padding + self.volume.shape[0], padding:padding +
            self.volume.shape[1], padding:padding + self.volume.shape[2]] =\
            self.volume
        self.volume = nest

        # lbls are all the indiviudal numbers that appear in the segmented file
        self.lbls = self.volume.copy()
        self.lbls = self.lbls.flatten()
        self.lbls = list(set(self.lbls))
        self.lbls.sort()
        self.lbls = self.lbls[1:] # get rid of 0 label

    def sequential_run(self):
        results = []
        q = Queue()
        get_loc_and_angle(self.volume, self.lbls, q, 0, self.padding,
                          self.debug, self.oblate)
    
        results.extend(q.get(True))
        
        results = np.array(results)
        np.save(self.fname + "_positions_orientation.npy", results)
        print("Found ", len(results), " particles")
    
    def processing(self):
        if self.force_sequential:
            return self.sequential_run()
        
        system_cores = mp.cpu_count()
        if self.specify_max_cores != 0:
            system_cores = min([self.specify_max_cores, system_cores])
        
        num_particles = len(self.lbls)
        chunk_length = num_particles // system_cores
        work_blocks = []
        
        for i in range(system_cores - 1):
            work_blocks.append(self.lbls[i*chunk_length: (i+1)*chunk_length])
        
        i = i + 1
        work_blocks.append(self.lbls[i*chunk_length:])
        
        
        workers = []
        results = []
        q = Queue()
        for i in range(system_cores):
            x =  Process(target=get_loc_and_angle,
                         args=(self.volume, work_blocks[i], q, i, self.padding,
                               self.debug, self.oblate))
            x.start()
            workers.append(x)
            
        for i in range(system_cores):
            results.extend(q.get(True))

        results = np.array(results)
        np.save(self.fname + "_positions_orientation.npy", results)
        print("Found ", len(results), " particles")
        
        for i in range(system_cores):
            workers[i].join()


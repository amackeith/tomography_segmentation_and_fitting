'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"
'''

from multiprocessing import Process, Queue
import multiprocessing as mp
import numpy as np
import scipy.ndimage as ndi


cpp_enabled = True
try:
    import watershedtools_cpp.watershedtools_cpp as watershedtools_cpp
except:
    print("CPP code not enabled")
    print("to build and enable it run ./build in the watershedtools_cpp folder")
    print("Proceeding using SLOWER python code")
    cpp_enabled = False


# in arr is a binary particle
def find_principal_axis(in_arr):
    arr = in_arr
    # com will be used as origin
    
    moment_of_inertia = np.zeros((3, 3)) * 1.0
    com = np.zeros(3) * 1.0
    
    if cpp_enabled:
        
        watershedtools_cpp.center_of_mass(arr, com)
        
        watershedtools_cpp.calculate_moment_of_inertia(in_arr,
                                                       moment_of_inertia,
                                                       com[0], com[1], com[2])

    else:
        com = np.array(ndi.measurements.center_of_mass(arr))
        shp = in_arr.shape
        for i in range(shp[0]):
            for j in range(shp[1]):
                for k in range(shp[2]):
                    if arr[i, j, k] == 1:
                        # x,y,z is the distance for the origin in each direction
                        x, y, z = np.array([i, j, k]) - com
                        moment_of_inertia[0, 0] = \
                            moment_of_inertia[0, 0] + y ** 2 + z ** 2  # I_xx
                        moment_of_inertia[1, 1] = \
                            moment_of_inertia[1, 1] + x ** 2 + z ** 2  # I_yy
                        moment_of_inertia[2, 2] = \
                            moment_of_inertia[2, 2] + x ** 2 + y ** 2  # I_zz
                        moment_of_inertia[0, 1] = \
                            moment_of_inertia[0, 1] - x * y  # I_xy
                        moment_of_inertia[1, 2] = \
                            moment_of_inertia[1, 2] - y * z  # I_yz
                        moment_of_inertia[0, 2] = \
                            moment_of_inertia[0, 2] - x * z  # I_xz
        
        moment_of_inertia[1, 0] = moment_of_inertia[0, 1]  # I_yx
        moment_of_inertia[2, 1] = moment_of_inertia[1, 2]  # I_zy
        moment_of_inertia[2, 0] = moment_of_inertia[0, 2]  # I_zx
        


    # such that eigen_values[i] is the ith Eval,
    # eigen_vectors[:,i] is the coresponding eigen vector
    eigen_values, eigen_vectors = np.linalg.eig(moment_of_inertia)
    
    # print "EVAL", eigen_values
    # print "EVEC\n",eigen_vectors
    
    # ident = np.array([[1.0,0,0], [0,1.0,0], [0,0,1.0]])

    # proof
    # print "should be zero", np.linalg.det(moment_of_inertia-ident*eigen_values[0]),np.linalg.det(moment_of_inertia-ident*eigen_values[1]),np.linalg.det(moment_of_inertia-ident*eigen_values[2])
    # print moment_of_inertia
    return eigen_values, eigen_vectors, com, moment_of_inertia


def get_phi_and_theta_from_principle_axis(pi_vec):
    # since theta in 0,pi and phi in 0, pi caputres every orientation
    # and since in those ranges r has a positive x component
    # if x dot r is < 0 use -r instead
    x_unit = np.array([1, 0, 0])
    z_unit = np.array([0, 0, 1])
    neg_y = np.array([0, -1, 0])
    
    if np.dot(x_unit, pi_vec) < 0:
        pi_vec = -1.0 * pi_vec

    # these projections use the euler angles to figure out what phi and theta are
    r_unit_vec = pi_vec / np.linalg.norm(pi_vec)
    theta = np.arccos(np.dot(z_unit, r_unit_vec))
    r_proj_on_xy = r_unit_vec.copy()
    r_proj_on_xy[2] = 0
    if np.linalg.norm(r_proj_on_xy) == 0.0:
        phi = 0.0
    else:
        r_proj_on_xy = r_proj_on_xy / np.linalg.norm(r_proj_on_xy)
        phi = np.arccos(np.dot(neg_y, r_proj_on_xy))
    
    return (phi, theta)



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
        
        val, vec, com, moment = find_principal_axis(arr)
        
        if oblate:  # particles are oblate
            index_of_max_eigen_val = np.argmax(val)
        else:  # prolate particle
            index_of_max_eigen_val = np.argmin(val)
        
        principle_axis = vec[:, index_of_max_eigen_val]
        phi, theta = get_phi_and_theta_from_principle_axis(principle_axis)
        
        correctd_com = com + np.array([x_avg, y_avg, z_avg]) - 2*padding
        # to correct for nesting it in the empty array
        
        # you need to take off 2 paddings, 1 for the padding of the whole volume
        # one for the padding of the subsection that is used when you calculate
        # the com and moment of inertia
        lentil = [correctd_com, np.array([phi, theta, 0]), principle_axis,
                  mask_vol, particle_index, (val, vec), moment]
        com_angle_list.append(lentil)
        
        if debug:
            print("core_num", core_num, "number", cnt,
                  "of", num, "mask vol", mask_vol)
        cnt = cnt + 1
    
    que.put(com_angle_list)


class orientations_from_moment_of_inertia:
    
    def __init__(self, fname, segmented_volume, padding=30,
                 oblate=True, debug=False, force_sequential=False,
                 specify_max_cores=0, threshold_of_binary=''):
        self.fname = fname[:-4]
        self.volume = np.int32(segmented_volume)
        self.padding = padding
        self.oblate = oblate
        self.debug = debug
        self.force_sequential = force_sequential
        self.specify_max_cores = specify_max_cores
        self.threshold_of_binary = threshold_of_binary
        
        # inset the volume with #padding voxels in each direction
        shp = np.array(self.volume.shape) + 2 * padding
        nest = np.zeros(shp)
        nest[padding:padding + self.volume.shape[0], padding:padding +
                                                             self.volume.shape[1],
        padding:padding + self.volume.shape[2]] = \
            self.volume
        self.volume = nest
        
        # lbls are all the indiviudal numbers that appear in the segmented file
        self.lbls = self.volume.copy()
        self.lbls = self.lbls.flatten()
        self.lbls = list(set(self.lbls))
        self.lbls.sort()
        self.lbls = self.lbls[1:]  # get rid of 0 label
    
    def sequential_run(self):
        results = []
        q = Queue()
        get_loc_and_angle(self.volume, self.lbls, q, 0, self.padding,
                          self.debug, self.oblate)
        
        results.extend(q.get(True))
        
        results = np.array(results)
        
        # see explination of positions_orientations.npy in processing()
        np.save(self.fname + self.threshold_of_binary +
                "_positions_orientation.npy", results)
        print("Found ", len(results), " particles")
    
    def processing(self):
        if self.force_sequential:
            return self.sequential_run()
        
        # since the granularity of the parallelization is so rough and there
        # is no load balancing I subtract one in the hope that no job is
        # deschedualled for very long even if there are other tasks being done
        # at the same time.
        system_cores = mp.cpu_count() - 1
        if self.specify_max_cores != 0:
            system_cores = min([self.specify_max_cores, system_cores])
        
        num_particles = len(self.lbls)
        chunk_length = num_particles // system_cores
        work_blocks = []
        
        for i in range(system_cores - 1):
            work_blocks.append(self.lbls[i * chunk_length: (i + 1) * chunk_length])
        
        i = i + 1
        work_blocks.append(self.lbls[i * chunk_length:])
        
        workers = []
        results = []
        q = Queue()
        for i in range(system_cores):
            x = Process(target=get_loc_and_angle,
                        args=(self.volume, work_blocks[i], q, i, self.padding,
                              self.debug, self.oblate))
            x.start()
            workers.append(x)
        
        for i in range(system_cores):
            results.extend(q.get(True))
        
        results = np.array(results, dtype=object)
        
        # this is to explain what each entry is. Every entry in this results
        # list is of the following form
        # particle = [center_of_mass, np.array([phi, theta, 0]),
        #             principle_axis, mask_vol, particle_index,
        #             (val, vec), moment]
        
        # center of mass is the center of mass with positions coords given by
        # the input volumes three axis.
        # phi and theta is the orientation of the principle axis (axis of
        # rotational symmetry of the given particle. To find it one uses this
        # rotation matrix:
        '''
        def make_rot_mat_simplified(phi, theta):
            rotmat = np.zeros((3, 3))
            rotmat[0, 0] = np.cos(phi)
            rotmat[0, 1] = np.sin(phi)
            rotmat[0, 2] = 0.
            
            rotmat[1, 0] = -np.cos(theta) * np.sin(phi)
            rotmat[1, 1] = np.cos(theta) * np.cos(phi)
            rotmat[1, 2] = np.sin(theta)
            
            rotmat[2, 0] = np.sin(theta) * np.sin(phi)
            rotmat[2, 1] = -np.sin(theta) * np.cos(phi)
            rotmat[2, 2] = np.cos(theta)
            
            return rotmat
        '''
        # on this kind of vector np.array([0,0,1.0]). The result of
        # principle_axis_from_phi_theta = np.dot(
        #                       make_rot_mat_simplified(phi, theta).T,
        #                       np.array([0,0,1.0])
        #### NOTE THE TRANSPOSE
        # should be the same within a sign of principle_axis.
        # the phi and theta are just there for making the particles as well
        # as showing the orientation with the two degrees of freedom it really has
        # mask_vol is the volume of the particular segment for this entry
        # particle_index is the segment label in the output volume from watershed pipeline
        # (val, vec) are the eigen values and vectors of the moment of inertia matrix
        # which is the final entry.
        
        np.save(self.fname + self.threshold_of_binary +
                "_positions_orientation.npy", results)
        print("Found ", len(results), " particles")
        
        for i in range(system_cores):
            workers[i].join()

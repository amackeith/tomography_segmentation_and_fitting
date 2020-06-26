'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"


This code parallelizes the lentil fitting process and is written of "oblate shapes" for prolate edit the line tagged with "#edit_here for prolate objects"
and chagne max to min
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


padding = 30 

def get_loc_and_angle(full_fitt, lbl_lst, que, core_num):
	full_fitt = full_fitt.copy()
	com_angle_list = []
	num = len(lbl_lst)
	cnt = 0

	for particle_index in lbl_lst:
		mask = np.array((full_fitt==particle_index)*1)
		mask_vol = np.sum(mask)
		np_where_info = np.where(mask==1)
		x_max = max(np_where_info[0])
		x_min = min(np_where_info[0])
		y_max = max(np_where_info[1])
		y_min = min(np_where_info[1])
		z_max = max(np_where_info[2])
		z_min = min(np_where_info[2])

		x_avg = (x_min+x_max)/2
		y_avg = (y_min+y_max)/2
		z_avg = (z_min+z_max)/2

		
		arr = mask[x_avg-padding:x_avg+padding, y_avg-padding:y_avg+padding, z_avg-padding:z_avg+padding]

		val,vec = find_principal_axis(arr)
		index_of_max_eigen_val = np.argmax(val) #edit_here for prolate objects
		principle_axis = vec[:,index_of_max_eigen_val]
		phi, theta = get_phi_and_theta_from_principle_axis(principle_axis)


		com = np.array(ndi.measurements.center_of_mass(mask))
		correctd_com = com-padding #to correct for nesting it in the empty array




		lentil = [correctd_com, np.array([phi,theta,0]), principle_axis, mask_vol, particle_index, (val, vec)]
		com_angle_list.append(lentil)

		

		print "core_num",core_num, "number", cnt, "of", num, "mask vol", mask_vol
		cnt = cnt+1

	que.put(com_angle_list)














if __name__ == '__main__':
	print "This is the grading copy"
	if len(sys.argv) != 3:
		print '\n\n'
		print sys.argv
		print "mom.py passed wrong number of arguments"
		exit()
	else:
		fname = sys.argv[1]
		title = sys.argv[2]


	data = np.load(fname)
	#tack on a buffer around the edges so you can slice a whole padding by padding by padding chunk out
	shp = np.array(data.shape)+2*padding
	nest = np.zeros(shp)
	nest[padding:padding+data.shape[0],padding:padding+data.shape[1],padding:padding+data.shape[2]] = data
	fitting = nest


	#lbls are all the indiviudal numbers that appear in the segmented file
	lbls = fitting.copy()
	lbls = lbls.flatten()
	lbls = list(set(lbls))
	

	lbls = map(float, lbls)
	lbls = map(int, lbls)
	lbls.sort()



	lbls = lbls[1:] #this gets rid of 0 label the background space

	start = time.time()

	num_particles = len(lbls)

	chunk_len = num_particles/6
	lbl1 = lbls[0:chunk_len]
	lbl2 = lbls[1*chunk_len:2*chunk_len]
	lbl3 = lbls[2*chunk_len:3*chunk_len]
	lbl4 = lbls[3*chunk_len:4*chunk_len]
	lbl5 = lbls[4*chunk_len:5*chunk_len]
	lbl6 = lbls[5*chunk_len:]


	q = Queue()
	results = []
	if mp.cpu_count()<7:
		print "using single core as there were less than 7 cores to use"
		get_loc_and_angle(fitting,lbls,q,0)

		results.extend(q.get(True))

	else:
		fitting1 = fitting.copy()
		fitting2 = fitting.copy()
		fitting3 = fitting.copy()
		fitting4 = fitting.copy()
		fitting5 = fitting.copy()
		fitting6 = fitting.copy()
		p1 = Process(target=get_loc_and_angle, args=(fitting1, lbl1, q, 1))
		p1.start()
		p2 = Process(target=get_loc_and_angle, args=(fitting2, lbl2, q,2))
		p2.start()
		p3 = Process(target=get_loc_and_angle, args=(fitting3, lbl3, q,3))
		p3.start()
		p4 = Process(target=get_loc_and_angle, args=(fitting4, lbl4, q,4))
		p4.start()
		
		p5 = Process(target=get_loc_and_angle, args=(fitting5, lbl5, q,5))
		p5.start()
		p6 = Process(target=get_loc_and_angle, args=(fitting6, lbl6, q,6))
		p6.start()
		
		for i in range(6):
			results.extend(q.get(True))



	results = np.array(results)
	np.save("safe_loc_norm_angle.npy", results)
	np.save("%s_loc_norm_angle.npy" % title, results)
	np.save("com_angle_list_%s_%s.npy" % (title,fname), results)
	print "Found %s lentils" % len(results)

	p1.join()
	p2.join()
	p3.join()
	p4.join()
	#p5.join()
	#p6.join()

	endTime = time.time()
	total = endTime - start
	
	print "The job took "+ str(total)+" seconds to complete"














import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import sys
from moment_of_inertia_orientation import find_principal_axis, get_phi_and_theta_from_principle_axis
import skimage.feature


#this is for visualization stuff
#from mpl_toolkits.mplot3d import Axes3D
#from mayavi import mlab



padding = 30 #this should be determined by the size of a lentil and how much padding it needs on each side of it


if len(sys.argv) != 2:
	print "mom.py passed wrong number of arguments"
	exit()
else:
	fname = sys.argv[1]


data = np.load(fname)
#tack on a buffer around the edges so you can slice a whole padding by padding by padding chunk out
shp = np.array(data.shape)+2*padding
nest = np.zeros(shp)
nest[padding:padding+data.shape[0],padding:padding+data.shape[1],padding:padding+data.shape[2]] = data
fitting = nest




#this means that nest[padding:-padding, padding:-padding, padding:-padding] is the actual fitting
#this means that x,y,z is com-padding

lbls = fitting.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()

print len(lbls)
com_angle_norm_vol = []
for particle_index in lbls[1:]:   #the first elt is the big part that is all zero
	print particle_index
	

	mask = np.array((fitting==particle_index)*1)
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
	index_of_max_eigen_val = np.argmax(val)
	principle_axis = vec[:,index_of_max_eigen_val]
	phi, theta = get_phi_and_theta_from_principle_axis(principle_axis)


	com = np.array(ndi.measurements.center_of_mass(mask))
	correctd_com = com-padding #to correct for nesting it in the empty array




	lentil = [correctd_com, np.array([phi,theta,0]), principle_axis, mask_vol]
	com_angle_norm_vol.append(lentil)


	#the rest is just for displaying stuff
	continue
	'''
	arr_com = np.array(ndi.measurements.center_of_mass(arr))
	tmp = original_binary_array[x_avg-padding:x_avg+padding, y_avg-padding:y_avg+padding, z_avg-padding:z_avg+padding]
	binary = tmp.copy()
	binary[binary>0]=1

	pi_scale = 10
	principle_axis = principle_axis*pi_scale
	
	#arr = np.complex(arr)

	#if mask_vol<1800 or (not np.array_equal(arr.shape, np.array([50,50,50]))):
		#continue

	#mlab.contour3d(arr)

	print lentil, particle_index
	print binary.shape, arr.shape

	mlab.pipeline.volume(mlab.pipeline.scalar_field(binary), color=(0,1,0))
	mlab.pipeline.volume(mlab.pipeline.scalar_field(arr), color=(1,0,0))
	mlab.quiver3d(arr_com[0],arr_com[1],arr_com[2],principle_axis[0],principle_axis[1],principle_axis[2], line_width=5, scale_factor=1,color=(0,0,1))
	mlab.show()
	'''
	

np.save("%s_com_angle_vol.npy" % fname, com_angle_norm_vol)
print "Found %s lentils" % len(com_angle_norm_vol)

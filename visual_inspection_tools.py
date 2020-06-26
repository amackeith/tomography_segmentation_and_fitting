#let this serve as a guide for visualiztion, just mess around with mayavi and see what you can do
#fname should always be a segmented file
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab

import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import sys
from moment_of_inertia_orientation import find_principal_axis, get_phi_and_theta_from_principle_axis
import skimage.feature



padding = 30 #this should be determined by the size of a lentil and how much padding it needs on each side of it

fname = "segmented_only_big_labels__file_in_lentil45fresh1210.npy_de_noised_de_ringed.npy_binary_thresh_0.62_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1000_.npy"

data = np.load(fname)
save_data = data.copy()
print data.shape
#tack on a buffer around the edges so you can slice a whole padding by padding by padding chunk out
shp = np.array(data.shape)+2*padding
nest = np.zeros(shp)
nest[padding:padding+data.shape[0],padding:padding+data.shape[1],padding:padding+data.shape[2]] = data
fitting = nest



#do the same to original array
#eroded binary array
original_binary_array = np.load("edm_of_ellipsiod_phase_file_in_lentil45fresh1210.npy_de_noised_de_ringed.npy_binary_thresh_0.62_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1000_.npy")
shp = np.array(original_binary_array.shape)+2*padding
nest2 = np.zeros(shp)
nest2[padding:padding+data.shape[0],padding:padding+data.shape[1],padding:padding+data.shape[2]] = original_binary_array
'''
com_etc = np.load("com_angle_list_small_slice.npy_segmented_only_big_labels__file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1500_.npy.npy")

com_etc = list(com_etc)
com_etc.sort(key=lambda x: x[4])



lbls = fitting.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()

vols = []


total_vol = []
other_vols = []
for i in lbls[1:]:
	i = int(i)
	print i, com_etc[i-1][4] 



	
	if 3400<com_etc[i-1][3]<30000:
		data[data==i]=0
		total_vol.append(com_etc[i-1][3])
	if com_etc[i-1][3]<3400:
		total_vol.append(com_etc[i-1][3])
		vols.append(com_etc[i-1][3])
		other_vols.append(com_etc[i-1][3])


print np.sum(vols), np.sum(total_vol), np.sum(vols)/(1.0*np.sum(total_vol))

plt.hist(other_vols, bins=100)
plt.show()

#mean is 2500
#this is for visualization stuff


def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=False):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    np.random.seed(100)



    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in xrange(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in xrange(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap




fitting = np.array(fitting>0)*1
nest2 = np.array(nest2>0)*1
arr = nest2-fitting

data_binary = np.array(data>0)*1
original_binary_array = np.array(original_binary_array>0)*1
arr = original_binary_array - data_binary

double_counts = save_data.copy()
double_counts[arr==0] = 0

# = 
np.save("double_counts.npy", double_counts)

colors = rand_cmap(1000)
for i in range(100):
	plt.title(str(i)+"helP")
	plt.imshow(double_counts[i], cmap=colors)
	plt.show()

from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab

arr = np.array(double_counts>0)*1

print arr.shape
mlab.pipeline.volume(mlab.pipeline.scalar_field(arr), color=(1,0,0))
mlab.show()


exit()
'''

lbls = fitting.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()

print len(lbls)
co=0
com_angle_norm_vol = []
for particle_index in lbls[1:]:   #the first elt is the big part that is all zero
	mask = np.array((fitting==particle_index)*1)
	mask_vol = np.sum(mask)

	if mask_vol>1000:
		continue
	else:
		co+=1
		print co
		continue
	

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

	#mlab.pipeline.volume(mlab.pipeline.scalar_field(binary), color=(0,1,0))
	mlab.pipeline.volume(mlab.pipeline.scalar_field(arr), color=(1,0,0))
	mlab.quiver3d(arr_com[0],arr_com[1],arr_com[2],principle_axis[0],principle_axis[1],principle_axis[2], line_width=5, scale_factor=5,color=(0,0,1))
	mlab.show()
	

print co

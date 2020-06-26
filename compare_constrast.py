import numpy as np
import matplotlib.pyplot as plt


arr = np.load("safe_loc_norm_angle.npy")
x = []
for i in arr:
    print x.append(i[3])
plt.hist(x, bins=200)
plt.show()
exit()


'''

tmp = np.load("edm_of_ellipsiod_phase_blur_file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.75_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy")





np_image_original = np.load("small_slice.npy_de_noised_de_ringed.npy")
rm_holes = np.load("rm_holes_file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy")
labeled_centers = np.load("labeled_centers_file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy")
segmented = np.load("segmented_using_grow_labels_file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy")
segmented_only_big_labels = np.load("segmented_only_big_labels__file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy")

'''
segmented = np.load("segmented_using_grow_labels_file_in_lentils45pour.npy_de_noised_de_ringed.npy_binary_thresh_0.67_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1500_.npy")
segmented_only_big_labels =np.load("segmented_using_grow_labels_file_in_lentils45pour.npy_de_noised_de_ringed.npy_binary_thresh_0.67_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1500_.npy")
print "ding"
lbls = segmented.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
#lbls.sort()
num_lentils_found1 = len(lbls)-1
print "ding"
lbls = segmented_only_big_labels.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
#lbls.sort()
num_lentils_found2 = len(lbls)-1
print "ding"
print num_lentils_found1, num_lentils_found2
exit()








segmented = np.array(segmented>0)*1
segmented_only_big_labels = np.array(segmented_only_big_labels>0)*1

print np.sum(segmented_only_big_labels - segmented)

for i in [470,]:

    plt.title(i)
    plt.imshow(segmented_only_big_labels[i], origin="lower")
    plt.show()





exit()



com_etc = np.load("com_angle_list_small_slice.npy_segmented_only_big_labels__file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1500_.npy.npy")
vols = []
total_vol = 0
bad_vol=0
print com_etc[44]
for i in com_etc[:]:
    print i[4]
    
    if i[3]>30000:
        print i[0]
        continue
    if i[3]>3430:
        total_vol += i[3]
        bad_vol = bad_vol+i[3]
    else:
        total_vol += i[3]
    vols.append(i[3])

print total_vol, bad_vol, bad_vol/(1.0*total_vol)

plt.hist(vols, bins=714)
plt.show()
exit()




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
            randRGBcolors[0] = [0, 0, 0]

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


lbls = segmented_only_big_labels.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()

max_vol = 4761
min_vol = 2846
'''
for particle_index in lbls[1:]:
    #print particle_index, len(lbls)
    mask = np.array(segmented_only_big_labels==particle_index)
    mask_vol = np.sum(mask)
    print particle_index, len(lbls)
    if not min_vol <mask_vol<max_vol:
        segmented_only_big_labels[mask] = 0
'''

#np.save("tmp.npy", segmented_only_big_labels)
segmented_only_big_labels_filt = np.load("tmp.npy")
color_map = rand_cmap(100)
save = False
length_of_movie = segmented_only_big_labels.shape[0]

segmented_only_big_labels = segmented_only_big_labels - segmented_only_big_labels_filt
for i in range(length_of_movie):
    plt.figure(figsize=(9,9))
    #if (i%10 != 0):
        #continue
    plt.imshow(segmented_only_big_labels[i], cmap=color_map)
    plt.show()



for i in range(length_of_movie):
	if (i%10 != 0):
		continue
	fig, ax = plt.subplots(3,3,figsize=(9,9))
	ax[0,0].imshow(np_image_original[i], cmap='gray', origin='lower')
	ax[0,0].set_title("Original slice %s \n(num found:%s) " % (i,len(lbls)))
	#x[0,1].imshow(thresh_np_image[i], cmap='gray', origin='lower')
	#ax[0,1].set_title("Binary Thresh %s" % threshold_of_binary)
	ax[0,2].imshow(rm_holes[i], cmap='gray', origin='lower')
	ax[0,2].set_title("Remove holes")
	#ax[1,0].imshow(eroded_binary[i], cmap='gray', origin='lower')
	#ax[1,0].set_title("Pre thresh on EDT")
	#ax[1,1].imshow(eroded_binary_thresh[i], cmap='gray', origin='lower')
	#ax[1,1].set_title("Post Thresh%s on EDT" % threshold_of_eroded_binary)
	ax[1,1].imshow(labeled_centers[i], origin='lower', cmap=color_map)
	ax[1,1].set_title("Labeled centers, \nThresh on EDT")

	#ax[1,2].imshow(edm_of_ellipsiod_phase_blur[i], cmap='gray', origin='lower')
	#ax[1,2].set_title("EDM,with guass filt sigma%s" % gauss_filt_sigma)
	#ax[2,0].imshow(edm_of_ellipsiod_phase[i], cmap='gray', origin='lower')
	#ax[2,0].set_title("EDM,thresh at %s" % threshold_of_edm2)
	
	ax[2,1].imshow(segmented[i], origin='lower', cmap=color_map)
	ax[2,1].set_title("Segments")
	ax[2,2].imshow(segmented_only_big_labels[i], origin='lower', cmap=color_map)
	ax[2,2].set_title("segmented only big labels vol")
	plt.tight_layout(h_pad=0.75)

	# plt.title(paramiters, y=1.08)
	#plt.savefig("freeze_frame_of"+paramiters+"_frame_%s.png" % i)
	plt.show()
	plt.clf()
	plt.close()


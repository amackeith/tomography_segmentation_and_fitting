import numpy as np 
import matplotlib.pyplot as plt


data = np.load("safe_loc_norm_angle.npy")
print data[33]

vols= []
for i in data[1:]:
    if i[3]<1000:
        vols.append(i[3])
print len(vols)
plt.hist(vols,bins=20)
plt.show()




exit()

arr = np.load("segmented_using_grow_labels_file_in_20181206lentils45_unworn.npy_de_noised_de_ringed.npy_binary_thresh_0.37_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1100_.npy")

ar = arr.copy()
ar = ar.flatten()
ar = list(set(ar))
vols = []
for i in ar[1:200]:
    print i
    s = np.sum(arr==i)
    vols.append(s)

plt.hist(vols, bins=40)
plt.show()





exit()

fname = "segmented_only_big_labels__file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy"
labeled = np.load(fname)
original_binary_array = np.load("rm_holes_file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy")
segmented = np.load("segmented_using_grow_labels_file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_1.4_min_vol_1500_.npy")


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



colors = rand_cmap(1000)
labeled = np.array(labeled>0)*1.0

segmented[labeled==1] = 0
original_binary_array = original_binary_array - labeled

for i in range(labeled.shape[0]):
	if i%3 ==0:
		print i

		plt.figure(figsize=(20,10))
		plt.title(i)
		#plt.imshow(labeled[i], cmap =colors, alpha=0.5)
		plt.imshow(segmented[i], cmap =colors, alpha=0.5)
		plt.imshow(original_binary_array[i], cmap='gray')
		#plt.imshow(labeled[i], cmap =colors, alpha=0.5)
		plt.imshow(segmented[i], cmap =colors, alpha=0.5)

		
		plt.show()
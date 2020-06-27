'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"


This is the meat of the process
'''

import matplotlib
import matplotlib.animation as animation
import random
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy
import scipy.ndimage as ndi
import sys
from hoshen_kopelmann_with_union_find import hoshen_kopelmann3d, grow_labels

if len(sys.argv) == 7:
    image_name_in = sys.argv[1]  # file name
    threshold_of_binary = float(sys.argv[2])  # space, lentil thresh
    threshold_of_eroded_binary = float(sys.argv[3])  # how far to cut away when looking for centers
    gauss_filt_sigma = float(sys.argv[4])  # smoothing
    threshold_of_edm2 = float(
        sys.argv[5])  # second threshold for determinin what is space what is not (and assigning height map)
    min_vol_for_segment = int(float(sys.argv[6]))  # anything smaller than this is thrown out
    print("label_watershed_method.py params:", image_name_in, threshold_of_binary, threshold_of_eroded_binary,
          gauss_filt_sigma, threshold_of_edm2, min_vol_for_segment)

else:
    print("wrong number of args", sys.argv)
    print(len(sys.argv))
    print(
        "should be python pgrm.py volume_file threshold_of_binary threshold_of_eroded_binary gauss_filt_sigma threshold_of_edm2 min_vol_for_segment")
    exit()

paramiters = "_file_in_%s_binary_thresh_%s_erosion_thresh_%s_gauss_filt_sigma_%s_edmthresh_%s_min_vol_%s_" % \
             (image_name_in, threshold_of_binary, threshold_of_eroded_binary, gauss_filt_sigma, threshold_of_edm2,
              min_vol_for_segment)

np.savetxt('paramiters.txt', [paramiters], fmt='%s')

np_image = np.load(image_name_in)
np_image_original = np_image[:]  # keep a copy of unmodified image for final comparison

# thresholding first time
thresh_np_image = np.array((np_image > threshold_of_binary) * 1.0)  # this threshold the image (take a look) with imshow
thresh_np_image = ndi.binary_erosion(thresh_np_image)
thresh_np_image = -1 * thresh_np_image + 1
np.save("thresh_np_image%s.npy" % paramiters, thresh_np_image)

# find out where the index of the biggest cluster (and set everything else to 1,
# this removes off voxels within the lentils)

hoshen_np_image = hoshen_kopelmann3d(thresh_np_image)
hoshen_np_image_flat = hoshen_np_image[:]
hoshen_np_image_flat = hoshen_np_image_flat.flatten()
counts = np.bincount(hoshen_np_image_flat)
background_lbl = np.argmax(counts)
rm_holes = np.int32(
    (hoshen_np_image != background_lbl) * 1)  # a threshold version with the holes filled in (important since we
# will be using eculedian distance maps)


np.save("hoshen_np_image%s.npy" % paramiters, hoshen_np_image)
np.save("rm_holes%s.npy" % paramiters, rm_holes)
hoshen_np_image = np.load("hoshen_np_image%s.npy" % paramiters)

rm_holes = np.load("rm_holes%s.npy" % paramiters)

# now calculate EDM and threshhold that map.
eroded_binary = ndi.morphology.distance_transform_edt(rm_holes)
eroded_binary_thresh = np.int8((eroded_binary > threshold_of_eroded_binary) * 1)
np.save("eroded_binary%s.npy" % paramiters, eroded_binary)
# using this threshheld map find the centers which should all be separated from eachother
labeled_centers = hoshen_kopelmann3d(eroded_binary_thresh)

np.save("labeled_centers%s.npy" % paramiters, labeled_centers)
labeled_centers = np.load("labeled_centers%s.npy" % paramiters)

# right branch, first take edm, then use gauss filter to get rid of non zero values in space
edm_of_ellipsiod_phase = ndi.morphology.distance_transform_edt(rm_holes)
edm_of_ellipsiod_phase = ndi.filters.gaussian_filter(edm_of_ellipsiod_phase, gauss_filt_sigma)
edm_of_ellipsiod_phase_blur = edm_of_ellipsiod_phase.copy()
np.save("edm_of_ellipsiod_phase_blur%s.npy" % paramiters, edm_of_ellipsiod_phase_blur)

edm_of_ellipsiod_phase[edm_of_ellipsiod_phase < threshold_of_edm2] = 0
np.save("edm_of_ellipsiod_phase%s.npy" % paramiters, edm_of_ellipsiod_phase)
edm_of_ellipsiod_phase = np.load("edm_of_ellipsiod_phase%s.npy" % paramiters)

print("starting grow labels")

# grow labels expands each label to take up the space that was previously threshheld, using the labeled centers
# found with hoshen kopelmann
segmented = grow_labels(edm_of_ellipsiod_phase, labeled_centers)
np.save("segmented_using_grow_labels%s.npy" % paramiters, segmented)

segmented = np.load("segmented_using_grow_labels%s.npy" % paramiters)

lbls = segmented.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()
num_lentils_found_seg = len(lbls) - 1

print("Found in first pass %s" % num_lentils_found_seg)

segmented_only_big_labels = segmented.copy()
# now remove those labels that are too small (less than half a lentil to be in there)


lbls = segmented_only_big_labels.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()

print("masking small labels")
cnt = 0
lbl_volume_list = []
volumes_list = []

for particle_index in lbls[1:]:
    mask = np.array(segmented_only_big_labels == particle_index)
    mask_vol = np.sum(mask)
    lbl_volume_list.append([particle_index, mask_vol])
    volumes_list.append(mask_vol)


####test
if False:
    plt.hist(volumes_list, bins=400)
    plt.show()
    min_vol_for_segment = int(float(input("Min Single Lent Value")))

min_vol_for_segment = 950.0

for i in lbl_volume_list:
    mask_vol = i[1]
    particle_index = i[0]
    mask = np.array(segmented_only_big_labels == particle_index)
    
    if mask_vol < min_vol_for_segment:
        segmented_only_big_labels[mask] = 0
        cnt = cnt + 1

print("finished masking small labels, starting grow labels for second time")
print("num lables removed", cnt)
edm_of_ellipsiod_phase = np.load("edm_of_ellipsiod_phase%s.npy" % paramiters)
segmented_only_big_labels = grow_labels(edm_of_ellipsiod_phase, segmented_only_big_labels)
np.save("segmented_only_big_labels_%s.npy" % paramiters, segmented_only_big_labels)

### vieng and recording


np_image_original = np.load(image_name_in)
segmented_only_big_labels = np.load("segmented_only_big_labels_%s.npy" % paramiters)
thresh_np_image = np.load("thresh_np_image%s.npy" % paramiters)
rm_holes = np.load("rm_holes%s.npy" % paramiters)
eroded_binary = np.load("eroded_binary%s.npy" % paramiters)
labeled_centers = np.load("labeled_centers%s.npy" % paramiters)
edm_of_ellipsiod_phase_blur = np.load("edm_of_ellipsiod_phase_blur%s.npy" % paramiters)
edm_of_ellipsiod_phase = np.load("edm_of_ellipsiod_phase%s.npy" % paramiters)
segmented = np.load("segmented_using_grow_labels%s.npy" % paramiters)


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
        print('Please choose "bright" or "soft" for type')
        return
    
    if verbose:
        print('Number of labels: ' + str(nlabels))
    
    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]
        
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
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]
        
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


colors = rand_cmap(1000, verbose=False)
seg_cmap = colors
# count how many different lentils there are

lbls = segmented_only_big_labels.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()
num_lentils_found = len(lbls) - 1

print("Found in only big %s" % num_lentils_found)

save = False
length_of_movie = segmented_only_big_labels.shape[0]

for i in range(length_of_movie):
    if (i % 10 != 0):
        continue
    fig, ax = plt.subplots(3, 3, figsize=(9, 9))
    ax[0, 0].imshow(np_image_original[i], cmap='gray', origin='lower')
    ax[0, 0].set_title("Original slice %s \n(num found:%s) " % (i, num_lentils_found))
    ax[0, 1].imshow(thresh_np_image[i], cmap='gray', origin='lower')
    ax[0, 1].set_title("Binary Thresh %s" % threshold_of_binary)
    ax[0, 2].imshow(rm_holes[i], cmap='gray', origin='lower')
    ax[0, 2].set_title("Remove holes")
    ax[1, 0].imshow(eroded_binary[i], cmap='gray', origin='lower')
    ax[1, 0].set_title("Pre thresh on EDT")
    # ax[1,1].imshow(eroded_binary_thresh[i], cmap='gray', origin='lower')
    # ax[1,1].set_title("Post Thresh%s on EDT" % threshold_of_eroded_binary)
    ax[1, 1].imshow(labeled_centers[i], origin='lower', cmap=seg_cmap)
    ax[1, 1].set_title("Labeled centers, \nThresh on EDT %s" % threshold_of_eroded_binary)
    
    ax[1, 2].imshow(edm_of_ellipsiod_phase_blur[i], cmap='gray', origin='lower')
    ax[1, 2].set_title("EDM,with guass filt sigma%s" % gauss_filt_sigma)
    ax[2, 0].imshow(edm_of_ellipsiod_phase[i], cmap='gray', origin='lower')
    ax[2, 0].set_title("EDM,thresh at %s" % threshold_of_edm2)
    
    ax[2, 1].imshow(segmented[i], origin='lower', cmap=seg_cmap)
    ax[2, 1].set_title("Segments")
    ax[2, 2].imshow(segmented_only_big_labels[i], origin='lower', cmap=seg_cmap)
    ax[2, 2].set_title("segmented only big labels vol %s" % min_vol_for_segment)
    plt.tight_layout(h_pad=0.75)
    
    # plt.title(paramiters, y=1.08)
    plt.savefig("freeze_frame_of" + paramiters + "_frame_%s.png" % i)
    plt.clf()
    plt.close()

if save == True:
    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title='Params:%s' % paramiters, artist='Matplotlib',
                    comment='Movie support!')
    writer = FFMpegWriter(fps=5, metadata=metadata)
    
    fig, ax = plt.subplots(3, 3, figsize=(9, 9))
    l, = plt.plot([], [], 'k-o')
    
    with writer.saving(fig, paramiters + ".mp4", length_of_movie):
        for i in range(length_of_movie):
            ax[0, 0].imshow(np_image_original[i], cmap='gray', origin='lower')
            ax[0, 0].set_title("Original slice %s \n(num found:%s) " % (i, num_lentils_found))
            ax[0, 1].imshow(thresh_np_image[i], cmap='gray', origin='lower')
            ax[0, 1].set_title("Binary Thresh %s" % threshold_of_binary)
            ax[0, 2].imshow(rm_holes[i], cmap='gray', origin='lower')
            ax[0, 2].set_title("Remove holes")
            ax[1, 0].imshow(eroded_binary[i], cmap='gray', origin='lower')
            ax[1, 0].set_title("Pre thresh on EDT")
            # ax[1,1].imshow(eroded_binary_thresh[i], cmap='gray', origin='lower')
            # ax[1,1].set_title("Post Thresh%s on EDT" % threshold_of_eroded_binary)
            ax[1, 1].imshow(labeled_centers[i], origin='lower', cmap=seg_cmap)
            ax[1, 1].set_title("Labeled centers, \nThresh on EDT %s" % threshold_of_eroded_binary)
            
            ax[1, 2].imshow(edm_of_ellipsiod_phase_blur[i], cmap='gray', origin='lower')
            ax[1, 2].set_title("EDM,with guass filt sigma%s" % gauss_filt_sigma)
            ax[2, 0].imshow(edm_of_ellipsiod_phase[i], cmap='gray', origin='lower')
            ax[2, 0].set_title("EDM,thresh at %s" % threshold_of_edm2)
            
            ax[2, 1].imshow(segmented[i], origin='lower', cmap=seg_cmap)
            ax[2, 1].set_title("Segments")
            ax[2, 2].imshow(segmented_only_big_labels[i], origin='lower', cmap=seg_cmap)
            ax[2, 2].set_title("segmented only big labels vol %s" % min_vol_for_segment)
            print(i, "of", length_of_movie, "Making label_watershed_method video")
            # plt.tight_layout(h_pad=0.75)
            # plt.title(paramiters, y=1.08)
            writer.grab_frame()
        # plt.clf()

print("Label Watershed Complete")

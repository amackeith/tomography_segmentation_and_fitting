'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"
'''

import sys, os, time
import numpy as np

time_start = time.time()

#Reconstructed tomography file name
fname = "shearworn45lentil20190110.recon.npy"
#if less than 2*width this will mask a cylinder centered at the selected centerx centery
mask_radius = 400
#binary thresh can be manually set here
threshold_of_binary = 0.69
#this is used to seperate the centers as "markers"
threshold_of_eroded_binary = 4.5
#this is just a blur that is used to smooth any odd discrepancies in the eculdian distance map
gauss_filt_sigma = 1.75
#this theresholds the blurred volume
threshold_of_edm2 = 2.0#1.4
#anything less than this volume is either incoperated into nearby stuff or thrown out
min_vol_for_segment = 1000


watershed_params = (threshold_of_binary, threshold_of_eroded_binary, gauss_filt_sigma, threshold_of_edm2, min_vol_for_segment)


#this is for the user to find the rings, usually they will be at 200, 200, but not always

if os.path.exists(fname):
	print "python denoise_lentils_22_select_center.py %s" % fname
	os.system("python denoise_lentils_22_select_center.py %s" % fname)
else:
	print "exit on find center"
	exit()




print "De Noising and De Ringing image"
print "python make_de_noised_de_ringed_given_cenxy.py %s" % fname
#given a centerx, centery (stored in a numpy file in the directory) this will denoise using Non local means and a tomography deringing function
os.system("python make_de_noised_de_ringed_given_cenxy.py %s %s" % (fname, mask_radius))


print "Use sliders to select threshold for binary"
if threshold_of_binary==None:
	os.system("python binary_threshold_selection.py %s_de_noised_de_ringed.npy" % fname)
	threshold_of_binary = np.load("%s_de_noised_de_ringed.npy_selected_binary_threshold.npy" % fname)
	threshold_of_binary = np.around(threshold_of_binary[0],2) #this is just to keep the file names reasonable


threshold_of_binary = np.around(threshold_of_binary,2)
#all instances of de_noised_de_ringed.npy must be replaced with "de_noised_de_ringed.npy" in order to implement denoising and derining

#print "de noising not done yet"

print "starting segmentation"
if os.path.exists("%s_de_noised_de_ringed.npy" % fname):
	###this line must be changed when de noising is integrated
	print "python label_watershed_method.py %s_de_noised_de_ringed.npy %s %s %s %s %s" % (fname, threshold_of_binary, threshold_of_eroded_binary, gauss_filt_sigma, threshold_of_edm2, min_vol_for_segment)
	os.system("python label_watershed_method.py %s_de_noised_de_ringed.npy %s %s %s %s %s" % (fname, threshold_of_binary, threshold_of_eroded_binary, gauss_filt_sigma, threshold_of_edm2, min_vol_for_segment))
else:
	print "De Noise failed"
	exit()


seg_name_params = "_file_in_%s_binary_thresh_%s_erosion_thresh_%s_gauss_filt_sigma_%s_edmthresh_%s_min_vol_%s_" % \
			(fname+"_de_noised_de_ringed.npy", threshold_of_binary, threshold_of_eroded_binary, gauss_filt_sigma, threshold_of_edm2, min_vol_for_segment)
seg_name = "segmented_only_big_labels_%s.npy" % seg_name_params


print "Finding centers and orientations"
if os.path.exists(seg_name):
	print "python moment_of_inertia_processing_multi_core.py %s %s" % (seg_name, fname)
	os.system("python moment_of_inertia_processing_multi_core.py %s %s" % (seg_name, fname))
else:
	print "watershed labeling failed"
	print seg_name
	#print "segmented_only_big_labels__file_in_small_slice.npy_de_noised_de_ringed.npy_binary_thresh_0.5_erosion_thresh_4.5_gauss_filt_sigma_1.75_edmthresh_2.0_min_vol_1500_.npy"
	exit()




#
#if os.system.exists("%s_com_angle_vol.npy" % seg_name) and visually_inspect:#
#	os.system("python ")


print "Mission Accomplished"

total = (time.time()-time_start)/60.0
print total, "minutes"
#starting visual inspection



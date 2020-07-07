'''
Author: Arthur Mackeith
Date: June 28, 2020


Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"

Expanded on by the Author with help from Kieran Murphy.
'''

import os, time
import tomography_preprocessing
import watershed_segmentation
from scipy.stats import scoreatpercentile



# parameters
# Reconstructed tomography file name
fname = "shearworn45lentil20190110.recon.npy"
outputfolder = os.path.expanduser("./trial_folder_cpp/")

# if less than 2*width this will mask a cylinder centered at the selected centerx centery
mask_radius = 400
# binary thresh can be manually set here
threshold_of_binary = 0.7401 #0.77 #None #0.69
# this is used to seperate the centers as "markers"
threshold_of_eroded_binary_percentile = 99.0 #4.5
# this is just a blur that is used to smooth any odd discrepancies in the eculdian distance map
gauss_filt_sigma = 1.75
# this theresholds the blurred volume
threshold_of_edm_percentile = 92.0 #2.0  # 1.4
# anything less than this volume is either incorporated into nearby stuff or thrown out
min_vol_for_segment = 1000
###
ring_artifact_in_center = True

debug = True

watershed_params = (threshold_of_binary,
                    threshold_of_eroded_binary_percentile,
                    gauss_filt_sigma,
                    threshold_of_edm_percentile,
                    min_vol_for_segment)


# takes some of the noise out of the image,
# removes ring artifact from tomography
# this is optimized for the equipment in the Jaeger Lab
# you may have to tweak parameters for other setups


def de_noise_and_de_ring_phase(output_folder, fname, ring_artifact_in_center, debug):
    if not os.path.exists(fname):
        raise FileNotFoundError(fname, " not found in directory")
    
    find_center = tomography_preprocessing.find_center(fname)
    
    # if it using the ring artifact is always in the center along the first
    # axis use
    # find_center.use_center_of_image()
    # if you need to select it use (you will need to make the artifact go
    # along the first axis to use this tool
    if ring_artifact_in_center:
        center_x, center_y = find_center.use_center_of_image()
    else:
        center_x, center_y = find_center.select_with_gui()
    
    preprocessor = tomography_preprocessing.tomo_preprocessor(
        outputfolder,
        fname,
        center_x,
        center_y,
        save_steps_and_images=debug,
        mask_radius=None)
    
    # using default for NL means kernel and ring width
    preprocessor.apply_non_local_means()
    preprocessor.apply_de_ring()
    
    # return the output array and (optionally save some output files)
    return preprocessor.save_output()

def watershed_phase(volume, fname, debug):
    binary_threshold = threshold_of_binary
    if threshold_of_binary is None:
        thresh_select = watershed_segmentation.threshold_selector(volume, fname, debug)
        binary_threshold = thresh_select.manualy_select_threshold()
    print("BEGIN")
    wsp = watershed_segmentation.watershed_pipeline(volume, binary_threshold,
                                              fname,
                                              threshold_of_eroded_binary_percentile,
                                              gauss_filt_sigma,
                                              threshold_of_edm_percentile,
                                              min_vol_for_segment,
                                              debug=True,
                                              use_better_labels=False)
    
    print(wsp.load_saved_parts())
    wsp.display_movie()
    exit()
    # wsp.remove_small_labels()
    # wsp.display_segments_step_through(1)
    # exit()
    ####
    s = time.time()
    print("one")
    wsp.binirize()
    print("two")
    wsp.remove_holes()
    print("three")
    wsp.find_centers_of_particles()
    print("four")
    wsp.create_grow_labels_edm()
    print("five")
    wsp.grow_labels()
    print("six")
    wsp.remove_small_labels()
    print("TOTAL TIME ON WATERSHED", (time.time() - s)/60, " Minutes")
    print("seven")
    wsp.display_standard_step_through()
    
def main():
    time_start = time.time()
    
    #de_ringed_volume = de_noise_and_de_ring_phase(outputfolder, fname, ring_artifact_in_center, debug)
    import numpy as np
    de_ringed_volume = \
        np.load("shearworn45lentil20190110.recon_preprocessed.npy")
    watershed_phase(de_ringed_volume,
                    outputfolder + fname, debug)
    
    print("Total time: ", (time.time() - time_start) / 60, " minutes")

if __name__ == "__main__":
    main()



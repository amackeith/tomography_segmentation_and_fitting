#!/usr/bin/env conda run -n seg_fit python
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
import particle_orientation
import argparse

# parameters
# Reconstructed tomography file name
fname = "shearworn45lentil20190110.recon.npy"
outputfolder = os.path.expanduser("./trial_folder_2/")


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


def watershed_phase(volume, fname, threshold_of_binary, gauss_filt_sigma,
                    threshold_of_eroded_binary_percentile,
                    threshold_of_edm_percentile, min_vol_for_segment,
                    debug):
    binary_threshold = threshold_of_binary
    if threshold_of_binary is None:
        thresh_select = \
            watershed_segmentation.threshold_selector(volume, fname, debug)
        binary_threshold = thresh_select.manualy_select_threshold()[0]
    print("BEGIN")
    wsp = watershed_segmentation.watershed_pipeline(volume, binary_threshold,
                                        fname,
                                        threshold_of_eroded_binary_percentile,
                                        gauss_filt_sigma,
                                        threshold_of_edm_percentile,
                                        min_vol_for_segment,
                                        debug=debug,
                                        use_better_labels=True)
    
    # wsp.load_saved_parts()
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
    print("TOTAL TIME ON WATERSHED", (time.time() - s) / 60, " Minutes")
    print("seven")
    # wsp.display_standard_step_through()
    wsp.display_movie()
    
    return wsp.segmented_only_big_labels


def center_of_mass_and_orientation_from_segment(fname, volume, padding=30,
                                                oblate=True, debug=True):
    part_orient = particle_orientation.orientations_from_moment_of_inertia(
        fname, volume, padding, oblate, debug)
    part_orient.processing()


def main():

    parser = argparse.ArgumentParser(description='Input for segmentation fit')
    
    parser.add_argument('-i', dest='fname', type=str,
                        help='npy input file (tomogram of packing)',
                        required=True)
    parser.add_argument('-o', dest='outputfolder', type=str,
                        default=str(int(time.time()*1000)),
                        help='output folder')
    parser.add_argument('-b', dest='threshold_of_binary',
                        default=None,
                        type=float,
                        help='threshold for eroded binary there is support to'
                             ' select this using a GUI if you do not include'
                             ' this flag')
    
    parser.add_argument('-g',
                        dest='gauss_filt_sigma', type=float,
                        default=1.75,
                        help='standard deviation of 3d gaussian filter used'
                             ' to blur the image slightly')
    
    parser.add_argument('-v',
                        dest='min_vol_for_segment',
                        default=1000,
                        type=int,
                        help='WARNING: this will be very dependent on shape : '
                             'Minimum size for a particle, all particles '
                             'less than this will be removed '
                             'and attempted to be re-assigned to other '
                             'particles')

    parser.add_argument('-r', dest='ring_artifact_in_center',
                        action='store_true',
                        help='The tomography ring artifact is located in the '
                             'center, if false you will need to select it on a '
                             'gui')
    parser.add_argument('-dont_save_files', dest='debug',
                        action='store_false',
                        help='If this flag is passed the only output will be '
                             'a few pngs and the final '
                             'orientation and center of'
                             ' mass list')
    parser.add_argument('-bp', dest='threshold_of_eroded_binary_percentile',
                        default='99.0',  type=float,
                        help=" this is used to seperate the "
                             "centers as 'markers'")
    parser.add_argument('-edm_p', dest='threshold_of_edm_percentile',
                        default='92.0', type=float,
                        help="this theresholds the blurred volume and sets "
                             "the voxels eligable to be part of the final"
                             " segments")
    
    args = parser.parse_args()
    print("WARNING: the default values will not work for every particle shape\n"
          "and tomography setup they must be adjusted to each project")
    threshold_of_binary = args.threshold_of_binary
    threshold_of_eroded_binary_percentile =\
        float(args.threshold_of_eroded_binary_percentile)
    threshold_of_edm_percentile =\
        float(args.threshold_of_edm_percentile)
    
    fname = args.fname
    gauss_filt_sigma = args.gauss_filt_sigma
    min_vol_for_segment = args.min_vol_for_segment
    outputfolder = os.path.expanduser(args.outputfolder)
    if outputfolder[-1] != "/":
        outputfolder += "/"
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
    
    ring_artifact_in_center = args.ring_artifact_in_center
    debug = args.debug
    
    watershed_params = (threshold_of_binary,
                        threshold_of_eroded_binary_percentile,
                        gauss_filt_sigma,
                        threshold_of_edm_percentile,
                        min_vol_for_segment)
    
    time_start = time.time()
    
    de_ringed_volume = de_noise_and_de_ring_phase(outputfolder,
     fname, ring_artifact_in_center, debug)
    import numpy as np
    #de_ringed_volume = \
    #    np.load("shearworn45lentil20190110.recon_preprocessed.npy")
    
    segmented_volume = watershed_phase(de_ringed_volume,
                                       outputfolder + fname,
                                       threshold_of_binary, gauss_filt_sigma,
                                       threshold_of_eroded_binary_percentile,
                                       threshold_of_edm_percentile,
                                       min_vol_for_segment,
                                       debug)
    
    center_of_mass_and_orientation_from_segment(fname, segmented_volume,
                                                padding=30, oblate=True,
                                                debug=debug)
    
    print("Total time: ", (time.time() - time_start) / 60, " minutes")


if __name__ == "__main__":
    main()

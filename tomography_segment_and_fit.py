#!/usr/bin/env python
'''
Author: Arthur Mackeith
Date: June 28, 2020


Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"

Expanded on by the Author with help from Kieran Murphy.


This code was re-written in July 2020 to run faster. It produces results
that are within 2.0e-8 voxels for position of particles and the orientation
vectors for the same particles have the properaty that:
np.linalg.norm(norm_original_code_orientation - this_code_orientation) < 1.0e-8. Which might be
attributable to changes between language versions or package versions.

'''

import numpy as np
import os, time
import tomography_preprocessing
import watershed_segmentation
import particle_orientation
import argparse



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
        output_folder,
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
    
    print("BEGIN watershed segmentation")
    wsp = watershed_segmentation.watershed_pipeline(volume, binary_threshold,
                                                    fname,
                                                    threshold_of_eroded_binary_percentile,
                                                    gauss_filt_sigma,
                                                    threshold_of_edm_percentile,
                                                    min_vol_for_segment,
                                                    debug=debug,
                                                    use_better_labels=True)
    
    # wsp.load_saved_parts()
    # return wsp.segmented_only_big_labels
    # wsp.display_segments_step_through(1)
    # exit()
    ####
    s = time.time()

    print("one: Binirize")
    wsp.binirize()
    print("two: Remove holes")
    wsp.remove_holes()
    print("three: find centers of particles")
    wsp.find_centers_of_particles()
    print("four: create Euclidean Distance Map to grow centers/labels")
    wsp.create_grow_labels_edm()
    print("five: Grow the centers/labels")
    wsp.grow_labels()
    print("six: Remove Small Labels")
    wsp.remove_small_labels()
    print("TOTAL TIME ON WATERSHED", (time.time() - s) / 60, " Minutes")
    print("seven: Create Display")
    # wsp.display_standard_step_through(10)
    wsp.display_movie()
    
    return wsp.segmented_only_big_labels


def center_of_mass_and_orientation_from_segment(fname, volume, padding=30,
                                                oblate=True, debug=True,
                                                force_sequential=False,
                                                threshold_of_binary=''):
    print("Calculating Center of Mass and Orientations from segments")
    part_orient = particle_orientation.orientations_from_moment_of_inertia(
        fname, volume, padding, oblate, debug,
        force_sequential=force_sequential,
        threshold_of_binary=threshold_of_binary)
    
    part_orient.processing()


def main():
    parser = argparse.ArgumentParser(description='Input for segmentation fit')
    
    parser.add_argument('-i', dest='fname', type=str,
                        help='npy input file (tomogram of packing)',
                        required=True)
    parser.add_argument('-o', dest='outputfolder', type=str,
                        default=str(int(time.time() * 1000)),
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
                        default='99.0', type=float,
                        help=" this is used to seperate the "
                             "centers as 'markers'")
    parser.add_argument('-edm_p', dest='threshold_of_edm_percentile',
                        default='92.0', type=float,
                        help="this theresholds the blurred volume and sets "
                             "the voxels eligable to be part of the final"
                             " segments")
    parser.add_argument('-prolate', dest='prolate',
                        action='store_true',
                        help="The default assumption is that the grains are"
                             " oblate this flag changes that to prolate, this "
                             "switches between the min and max inertial axis"
                             " being used as the orientation vector.")
    
    parser.add_argument('-de_noised_volume', dest='de_noised_volume',
                        action='store_true',
                        help="skip de noising and use fname this volume instead")
    
    args = parser.parse_args()
    print("WARNING: the default values will not work for every particle shape\n"
          "and tomography setup they must be adjusted to each project")
    threshold_of_binary = args.threshold_of_binary
    threshold_of_eroded_binary_percentile = \
        float(args.threshold_of_eroded_binary_percentile)
    threshold_of_edm_percentile = \
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
    oblate = not args.prolate
    de_noised_volume = args.de_noised_volume

    watershed_params = (threshold_of_binary,
                        threshold_of_eroded_binary_percentile,
                        gauss_filt_sigma,
                        threshold_of_edm_percentile,
                        min_vol_for_segment)
    
    time_start = time.time()
    
    if de_noised_volume == True:
        print("SKIPPING DE NOSING")
        de_ringed_volume = np.load(fname)
    else:
        de_ringed_volume = de_noise_and_de_ring_phase(outputfolder,
                                                      fname, ring_artifact_in_center, debug)
    
    
    segmented_volume = watershed_phase(de_ringed_volume,
                                       outputfolder + fname,
                                       threshold_of_binary, gauss_filt_sigma,
                                       threshold_of_eroded_binary_percentile,
                                       threshold_of_edm_percentile,
                                       min_vol_for_segment,
                                       debug)
    
    center_of_mass_and_orientation_from_segment(
        outputfolder + fname,
        segmented_volume,
        padding=30, oblate=oblate,
        debug=debug,
        threshold_of_binary=str(threshold_of_binary))

    print("Total time: ", (time.time() - time_start) / 60, " minutes")


if __name__ == "__main__":
    main()

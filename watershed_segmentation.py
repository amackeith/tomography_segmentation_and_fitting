'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"


Gui to help select threshold,
basically this is just three sliders using matplot lib, an example to match
and it saves the threshold to 2 decimal places
'''


from matplotlib.widgets import Slider
import matplotlib.animation as animation
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from scipy.stats import scoreatpercentile
from watershedtools_py import \
    hoshen_kopelmann3d, grow_labels, better_labels

# this block is to allow the use of this code without cpp speedup
# note this code reduces the runtime of several steps by orders of magnitude
# but you only need to run the code a few times and don't have cpp configured
# and don't know how, the python code is just as correct.
cpp_enabled = True
try:
    import watershedtools_cpp.watershedtools_cpp as watershedtools_cpp
except:
    print("CPP code not enabled")
    print("to build and enable it run ./build in the watershedtools_cpp folder")
    print("Proceeding using SLOWER python code")
    cpp_enabled = False


def rand_cmap(nlabels, kind='bright',
                  first_color_black=True,
                  last_color_black=False, verbose=False):
        """
        Creates a random colormap to be used together with
         matplotlib. Useful for segmentation tasks
        :param nlabels: Number of labels (size of colormap)
        :param kind: 'bright' for strong colors, 'soft' for pastel colors
        :param first_color_black: Option to use first color as black, True or False
        :param last_color_black: Option to use last color as black, True or False
        :param verbose: Prints the number of labels and shows the colormap. True or False
        :return: colormap for matplotlib
        """

        # this function from  https://github.com/delestro/rand_cmap
        from matplotlib.colors import LinearSegmentedColormap
        import colorsys
        import numpy as np

        np.random.seed(100)

        if kind not in ('bright', 'soft'):
            print('Please choose "bright" or "soft" for kind')
            return

        if verbose:
            print('Number of labels: ' + str(nlabels))

        # Generate color map for bright colors, based on hsv
        if kind == 'bright':
            randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                              np.random.uniform(low=0.2, high=1),
                              np.random.uniform(low=0.9, high=1))
                             for i in range(nlabels)]

            # Convert HSV list to RGB
            randRGBcolors = []
            for HSVcolor in randHSVcolors:
                randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0],
                                                         HSVcolor[1],
                                                         HSVcolor[2]))

            if first_color_black:
                randRGBcolors[0] = [0, 0, 0]

            if last_color_black:
                randRGBcolors[-1] = [0, 0, 0]

            random_colormap = LinearSegmentedColormap.from_list('new_map',
                                                                randRGBcolors,
                                                                N=nlabels)

        # Generate soft pastel colors, by limiting the RGB spectrum
        if kind == 'soft':
            low = 0.6
            high = 0.95
            randRGBcolors = [(np.random.uniform(low=low, high=high),
                              np.random.uniform(low=low, high=high),
                              np.random.uniform(low=low, high=high))
                             for i in range(nlabels)]

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

            cb = colorbar.ColorbarBase(ax, cmap=random_colormap,
                                       norm=norm,
                                       spacing='proportional', ticks=None,
                                       boundaries=bounds,
                                       format='%1i', orientation=u'horizontal')

        return random_colormap


class threshold_selector:
    
    def __init__(self, input_array, fname=None, debug=False):
        self.de_noised_volume = input_array
        self.fname = fname
        self.debug = debug
    
    def manualy_select_threshold(self):
        self.histvals = self.de_noised_volume.copy()
        self.histvals = self.histvals.flatten()
        
        self.fig, self.ax = plt.subplots(1, 3, figsize=(16, 8))
        self.ax[0].hist(self.histvals, bins=100)
        (n, bins, patches) = plt.hist(self.histvals, bins=100)
        ymax = np.max(n[10:]) * 1.05
        self.ax[0].set_xlim(xmin=0.05)
        self.ax[0].set_ylim(ymax=ymax)
        self.ax[0].set_title("Tomography value histogram (not including 0)")
        self.fig.subplots_adjust(left=0.1, bottom=0.25)
        min0 = 0
        max0 = (self.de_noised_volume.shape[0] - 1) // 2
        im1 = self.ax[1].imshow(self.de_noised_volume[max0])
        self.ax[1].set_title("Select Threshold using slider")
        
        example = np.load("example_binary_thresh.npy")
        min1 = 0
        max1 = (example.shape[0] - 1) // 2
        im2 = self.ax[2].imshow(example[max1])
        
        axcolor = 'lightgoldenrodyellow'
        axmin = self.fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
        axmax = self.fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        ax_example = self.fig.add_axes([0.25, 0.05, 0.65, 0.03],
                                       facecolor=axcolor)
        
        half = (np.min(bins) + np.max(bins)) // 2
        threshold_slider = Slider(axmin, 'Threshold',
                                  np.min(bins), np.max(bins), valinit=half)
        slice_slider = Slider(axmax,
                              'Slice', 0,
                              self.de_noised_volume.shape[0] - 1, valinit=max0)
        example_slider = Slider(ax_example, 'Example Slice',
                                0, example.shape[0] - 1, valinit=max1)
        
        def update(val):
            self.ax[1].cla()
            self.ax[1].imshow(
                np.array(
                    (self.de_noised_volume > threshold_slider.val)
                    * 1.0)[int(slice_slider.val)])
            self.ax[2].cla()
            self.ax[2].imshow(example[int(example_slider.val)])
            self.fig.canvas.draw()
        
        threshold_slider.on_changed(update)
        slice_slider.on_changed(update)
        example_slider.on_changed(update)
        plt.show()
        
        thresh = np.array([threshold_slider.val])
        
        if self.debug:
            np.save("%s_selected_binary_threshold.npy" % self.fname[:-4],
                    thresh)
        
        return thresh

        # TODO: implement auto threshold selection if tomography's and denoising
        # TODO: becomes reliable enough
        # TODO: e.g. have the threshold be on the intersect of two gaussians
        # TODO: fit to the density distribution of the volume (not includign 0)


class watershed_pipeline:
    # for the intended use of this the stages are called in order that they
    # are defined in this class
    
    def __init__(self, de_noised_volume, threshold_of_binary,
                 fname,
                 threshold_of_eroded_binary_percentile,
                 gauss_filt_sigma,
                 threshold_of_edm_percentile,
                 min_vol_for_segment,
                 debug=False,
                 use_better_labels=False):
        
        # parameters
        self.de_noised_volume = de_noised_volume
        self.threshold_of_binary = threshold_of_binary
        self.threshold_of_eroded_binary_percentile = threshold_of_eroded_binary_percentile
        self.fname = fname[:-4] + str(threshold_of_binary)
        self.gauss_filt_sigma = gauss_filt_sigma
        self.threshold_of_edm_percentile = threshold_of_edm_percentile
        self.min_vol_for_segment = min_vol_for_segment
        self.debug = debug
        self.use_better_labels = use_better_labels
        
        # keep a copy of unmodified image for final comparison
        if self.de_noised_volume is not None:
            self.original_volume = self.de_noised_volume.copy()

        self.parameters = "%s %s %s %s %s %s" % \
                     (fname, threshold_of_binary,
                      threshold_of_eroded_binary_percentile, gauss_filt_sigma,
                      threshold_of_edm_percentile,
                      min_vol_for_segment)
        if self.debug:
            
    
            np.savetxt(self.fname + '_watershed_parameters.txt',
                       [self.parameters], fmt='%s')
            
            np.save(self.fname + "_de_noised.npy", self.original_volume)

        
        
        
        
    
        
        
    def binirize(self):
        # thresholding first time
        self.thresh_np_image = np.array(
            (self.de_noised_volume > self.threshold_of_binary) * 1.0)
        # this threshold the image (take a look) with imshow
        self.thresh_np_image = ndi.binary_erosion(self.thresh_np_image)
        self.thresh_np_image = -1 * self.thresh_np_image + 1
        self.thresh_np_image = np.int32(self.thresh_np_image)
        if self.debug:
            np.save(self.fname + "_binary.npy", self.thresh_np_image)
        
    def remove_holes(self):
        # remove any holes inside so
        # erosion doesn't make holes appear inside the
        # particles
        if cpp_enabled:
            
            tmp = np.int32(self.thresh_np_image != 0)
            watershedtools_cpp.hoshen_kopelman3d_interface(tmp, self.debug)
            self.hoshen_np_image = tmp

        else:
            self.hoshen_np_image = hoshen_kopelmann3d(self.thresh_np_image)
            
        self.hoshen_np_image_flat = self.hoshen_np_image[:]
        self.hoshen_np_image_flat = self.hoshen_np_image_flat.flatten()
        counts = np.bincount(self.hoshen_np_image_flat)
        background_lbl = np.argmax(counts)
        self.rm_holes = np.int32(
            (self.hoshen_np_image != background_lbl) * 1)
        # a threshold version with the holes filled in (important since we
        # will be using eculedian distance maps)
        if self.debug:
            np.save(self.fname + "_rm_holes.npy", self.rm_holes)

    def find_centers_of_particles(self):
        # calculate Euclidian Distane map and threshhold that map to
        # isolate the centers of grains.
        self.eroded_binary =\
            ndi.morphology.distance_transform_edt(self.rm_holes)
        
        min_p, max_p = \
            scoreatpercentile(self.eroded_binary,
                              (self.threshold_of_eroded_binary_percentile, 99.5))

        self.eroded_binary_thresh = np.int8((self.eroded_binary > min_p) * 1)
        # using this threshheld map find the
        # centers which should all be separated from eachother
        
        # Hoshen Kopelmann then actually labels all of the isolated segments
        if cpp_enabled:
            self.labeled_centers = np.int32(self.eroded_binary_thresh.copy())
            watershedtools_cpp.hoshen_kopelman3d_interface(self.labeled_centers,
                                                           self.debug)
        else:
            self.labeled_centers = hoshen_kopelmann3d(self.eroded_binary_thresh)
            
        if self.use_better_labels:
            self.labeled_centers = better_labels(self.labeled_centers)

        if self.debug:
            np.save(self.fname + "_eroded_binary.npy", self.eroded_binary)
            
            
            np.save(self.fname + "_labeled_centers.npy",
                    self.labeled_centers)
            
    def create_grow_labels_edm(self):
        # this part creates the "topology" of the "watershed"
        # To use the labels found in find centesr of
        # particles, we must grow them
        # until they include the entire (or almost entire) particle
        # to do this we need a way to grow them, this is achieved by taking
        # a euclidian distance map of the binary packing again, (and in a later
        # step) expanding the labels along the positive gradient of the EDM.

        # right branch, first take edm, then use gauss
        # filter to get rid of non zero values in space
        self.edm_of_ellipsiod_phase = \
            ndi.morphology.distance_transform_edt(self.rm_holes)
        
        # this blur is to make the gradients smoother (tho I have considerd
        # removign it.
        self.edm_of_ellipsiod_phase = \
            ndi.filters.gaussian_filter(self.edm_of_ellipsiod_phase,
                                        self.gauss_filt_sigma)
        self.edm_of_ellipsiod_phase_blur = self.edm_of_ellipsiod_phase.copy()
        
        if self.debug:
            np.save(self.fname + "_edm_blur.npy",
                    self.edm_of_ellipsiod_phase_blur)

        if self.threshold_of_edm_percentile != -1:
            min_p, max_p = \
                scoreatpercentile(self.edm_of_ellipsiod_phase,
                                  (self.threshold_of_edm_percentile, 99.5))
            self.edm_of_ellipsiod_phase[
                        self.edm_of_ellipsiod_phase < min_p] = 0

            self.edm_of_ellipsiod_phase[
                self.edm_of_ellipsiod_phase < 2.0] = 0
        else:
            # In this case the -1 is a flag which means we want to grow to the 
            # original size of the binirized volume
            self.edm_of_ellipsiod_phase[self.rm_holes == 0] = 0
        # HEREERERERER
        
        

        if self.debug:
            np.save(self.fname + "_edm_blur_threshold.npy",
                    self.edm_of_ellipsiod_phase)


    def grow_labels(self):
        # this usually results in fractured cells, so you must merge the smaller
        # ones
        if cpp_enabled:
            # done array is an array used inside grow_labels_initerface
            # I just initiliaze it out here.
            done_array = np.zeros(self.edm_of_ellipsiod_phase.shape)
            self.segmented = self.labeled_centers.copy()*1.0
            watershedtools_cpp.grow_labels_interface(done_array,
                                                 self.segmented,
                                                 self.edm_of_ellipsiod_phase,
                                                 self.debug)
            
            
        else:
            self.segmented = grow_labels(self.edm_of_ellipsiod_phase,
                                         self.labeled_centers)
        
        if self.use_better_labels:
            self.segmented = better_labels(self.segmented)
        
        if self.debug:
            np.save(self.fname + "_first_grow_labels.npy", self.segmented)
            
    def remove_small_labels(self):
        self.segmented_only_big_labels = self.segmented.copy()
        # now remove those labels that are too small (less than half
        # a lentil to be in there)
        flt = self.segmented_only_big_labels.copy()
        flt = flt.flatten()
        size_dict = {}
        
        for label in flt:
            if label != 0:
                if label in size_dict:
                    size_dict[label] = size_dict[label] + 1
                else:
                    size_dict[label] = 1
    
    
        lbls = self.segmented_only_big_labels.copy()
        lbls = lbls.flatten()
        lbls = list(set(lbls))
        lbls.sort()
        self.num_lentils_found_seg = len(lbls) - 1
        print("Found in first pass %s" % self.num_lentils_found_seg)

        cnt = 0

        lbl_volume_list = []
        volumes_list = []
        
        for lbl in size_dict.keys():
            mask_vol = size_dict[lbl]
            lbl_volume_list.append([lbl, mask_vol])
            volumes_list.append(mask_vol)

        if self.debug:
            plt.clf()
            plt.hist(volumes_list, bins=400)
            plt.savefig(self.fname + "_label_volume_distribution.png")
            


        if cpp_enabled:
            self.segmented_only_big_labels = np.int32(self.segmented_only_big_labels)
            watershedtools_cpp.remove_small_labels_interface(
                self.min_vol_for_segment,
                self.segmented_only_big_labels)

        else:
                
            for i in lbl_volume_list:
                mask_vol = i[1]
                particle_index = i[0]
                mask = np.array(self.segmented_only_big_labels == particle_index)
        
                if mask_vol < self.min_vol_for_segment:
                    self.segmented_only_big_labels[mask] = 0
                    cnt = cnt + 1

            
            print("num lables removed", cnt)

        print("finished masking small labels, "
              "starting grow labels for second time")
        if cpp_enabled:
            # done array is an array used inside grow_labels_initerface
            # I just initiliaze it out here.
            done_array = np.zeros(self.edm_of_ellipsiod_phase.shape)
            self.segmented_only_big_labels = self.segmented_only_big_labels*1.0
            watershedtools_cpp.grow_labels_interface(done_array,
                                                 self.segmented_only_big_labels,
                                                 self.edm_of_ellipsiod_phase,
                                                 self.debug)
        else:
            self.segmented_only_big_labels = \
                grow_labels(self.edm_of_ellipsiod_phase,
                            self.segmented_only_big_labels)
            
        lbls = self.segmented_only_big_labels.copy()
        lbls = lbls.flatten()
        lbls = list(set(lbls))
        lbls.sort()
        self.num_lentils_found = len(lbls) - 1
        print('Second pass, number particles found', self.num_lentils_found)
        
        if self.use_better_labels:
            self.segmented_only_big_labels = \
                better_labels(self.segmented_only_big_labels)
        
        #if self.debug:

        np.save(self.fname + "_second_grow_labels.npy",
                    self.segmented_only_big_labels)
    
    
    def load_saved_parts(self):
        cnt = 0
        try:
            try: #the first one could be in two different places
                self.original_volume = np.load(self.fname + "_preprocessed.npy")
                cnt += 1
            except FileNotFoundError:
                try:
                    self.original_volume = \
                        np.load(self.fname + "_de_noised.npy")
                    cnt += 1
                except FileNotFoundError:
                    nm = self.fname[:-len(str(self.threshold_of_binary))]
                    self.original_volume = \
                        np.load( nm + "_original.npy")
                    cnt += 1

            try:
                self.thresh_np_image = np.load(self.fname + "_binary.npy")
                cnt += 1
            except FileNotFoundError as e:
                print(e)
            try:
                self.rm_holes = np.load(self.fname + "_rm_holes.npy")
                cnt += 1
            except FileNotFoundError as e:
                print(e)
            try:
                self.eroded_binary = np.load(self.fname + "_eroded_binary.npy")
                cnt += 1
            except FileNotFoundError as e:
                print(e)
            try:
                self.labeled_centers = \
                    np.load(self.fname + "_labeled_centers.npy")
                cnt += 1
            except FileNotFoundError as e:
                print(e)
            try:
                self.edm_of_ellipsiod_phase_blur = \
                    np.load(self.fname + "_edm_blur.npy")
                cnt += 1
            except FileNotFoundError as e:
                print(e)
            try:
                self.edm_of_ellipsiod_phase = \
                    np.load(self.fname + "_edm_blur_threshold.npy")
                cnt += 1
            except FileNotFoundError as e:
                print(e)
            try:
                self.segmented = np.load(self.fname + "_first_grow_labels.npy")
                cnt += 1
            except FileNotFoundError as e:
                print(e)
            try:
                self.segmented_only_big_labels = \
                    np.load(self.fname + "_second_grow_labels.npy")
                lbls = self.segmented_only_big_labels.copy()
                lbls = lbls.flatten()
                lbls = list(set(lbls))
                self.num_lentils_found = len(lbls) - 1

                cnt += 1
            except FileNotFoundError as e:
                print(e)

        except FileNotFoundError as e:
            print(e)
        finally:
            return cnt == 9, cnt
        
    def display_standard_step_through(self, step_size=1):
        
        #this is to make the color map consistent accross slices
        ones_with_color = [self.labeled_centers, self.segmented, self.segmented_only_big_labels]
        for c in ones_with_color:
            c[:, 0, 0] = 10000
            
        colors = rand_cmap(10000, verbose=False)
        seg_cmap = colors
        
        ## this section (at the cost of 1 bad px makes the color
        # maps consistent from shot to shot
        # so this step should only be done after saving
        self.labeled_centers[:,0,0] = 10000
        self.segmented[:, 0, 0] = 10000
        self.segmented_only_big_labels[:, 0, 0] = 10000

        length_of_movie = self.segmented_only_big_labels.shape[0]
        
        for i in range(length_of_movie):
            if (i % step_size != 0):
                continue
            fig, ax = plt.subplots(3, 3, figsize=(9, 9))
            ax[0, 0].imshow(self.original_volume[i],
                            cmap='gray', origin='lower')
            ax[0, 0].set_title("Original slice %s \n(num found:%s) "
                               % (i, self.num_lentils_found))
            ax[0, 1].imshow(self.thresh_np_image[i],
                            cmap='gray', origin='lower')
            ax[0, 1].set_title("Binary Thresh %s"
                               % self.threshold_of_binary)
            ax[0, 2].imshow(self.rm_holes[i], cmap='gray', origin='lower')
            ax[0, 2].set_title("Remove holes")
            ax[1, 0].imshow(self.eroded_binary[i], cmap='gray', origin='lower')
            ax[1, 0].set_title("Pre thresh on EDT")
            ax[1, 1].imshow(self.labeled_centers[i],
                            origin='lower', cmap=seg_cmap, interpolation='none')
            ax[1, 1].set_title("Labeled centers, \nThresh on EDT %s"
                               % self.threshold_of_eroded_binary_percentile)
    
            ax[1, 2].imshow(self.edm_of_ellipsiod_phase_blur[i],
                            cmap='gray', origin='lower')
            ax[1, 2].set_title("EDM,with guass filt sigma%s"
                               % self.gauss_filt_sigma)
            ax[2, 0].imshow(self.edm_of_ellipsiod_phase[i],
                            cmap='gray', origin='lower')
            ax[2, 0].set_title("EDM,thresh at %s" %
                               self.threshold_of_edm_percentile)
    
            ax[2, 1].imshow(self.segmented[i],
                            origin='lower', cmap=seg_cmap, interpolation='none')
            ax[2, 1].set_title("Segments")
            ax[2, 2].imshow(self.segmented_only_big_labels[i],
                            origin='lower', cmap=seg_cmap, interpolation='none')
            ax[2, 2].set_title("segmented only big labels vol %s"
                               % self.min_vol_for_segment)
            plt.tight_layout(h_pad=0.75)
    
            # plt.title(paramiters, y=1.08)
            plt.savefig(self.fname + "_frame_%s.png" % i)
            plt.clf()
            plt.close()

        for c in ones_with_color:
            c[:, 0, 0] = 0

    def display_movie(self):
    
        # this is to make the color map consistent accross slices
        ones_with_color = [self.labeled_centers, self.segmented, self.segmented_only_big_labels]
        for c in ones_with_color:
            c[:, 0, 0] = 10000
        
        colors = rand_cmap(10000, verbose=False)
        seg_cmap = colors
        FFMpegWriter = animation.writers['ffmpeg']
        
        metadata = dict(title='Params:%s' % self.parameters, artist='Matplotlib',
                        comment='Movie support!')
        writer = FFMpegWriter(fps=5, metadata=metadata)
    
        fig, ax = plt.subplots(3, 3, figsize=(9, 9))
        l, = plt.plot([], [], 'k-o')
        
        length_of_movie = self.original_volume.shape[0]
        with writer.saving(fig, self.fname + ".mp4", length_of_movie):
            for i in range(length_of_movie):
                
                ###### without this block of code this process will take like
                ###### 45 minutes
                for ii in range(3):
                    for jj in range(3):
                        ax[ii,jj].clear()
                ######
                        
                ax[0, 0].imshow(self.original_volume[i],
                                cmap='gray', origin='lower')
                ax[0, 0].set_title("Original slice %s \n(num found:%s) "
                                   % (i, self.num_lentils_found))
                ax[0, 1].imshow(self.thresh_np_image[i],
                                cmap='gray', origin='lower')
                ax[0, 1].set_title("Binary Thresh %s"
                                   % self.threshold_of_binary)
                ax[0, 2].imshow(self.rm_holes[i],
                                cmap='gray', origin='lower')
                ax[0, 2].set_title("Remove holes")
                ax[1, 0].imshow(self.eroded_binary[i],
                                cmap='gray', origin='lower')
                ax[1, 0].set_title("Pre thresh on EDT")
                ax[1, 1].imshow(self.labeled_centers[i],
                                origin='lower',
                                cmap=seg_cmap, interpolation='none')
                ax[1, 1].set_title("Labeled centers, \nThresh on EDT %s"
                                   % self.threshold_of_eroded_binary_percentile)
    
                ax[1, 2].imshow(self.edm_of_ellipsiod_phase_blur[i],
                                cmap='gray', origin='lower')
                ax[1, 2].set_title("EDM,with guass filt sigma%s"
                                   % self.gauss_filt_sigma)
                ax[2, 0].imshow(self.edm_of_ellipsiod_phase[i],
                                cmap='gray', origin='lower')
                ax[2, 0].set_title("EDM,thresh at %s" %
                                   self.threshold_of_edm_percentile)
    
                ax[2, 1].imshow(self.segmented[i],
                                origin='lower',
                                cmap=seg_cmap, interpolation='none')
                ax[2, 1].set_title("Segments")
                ax[2, 2].imshow(self.segmented_only_big_labels[i],
                                origin='lower',
                                cmap=seg_cmap, interpolation='none')
                ax[2, 2].set_title("segmented only big labels vol %s"
                                   % self.min_vol_for_segment)
                #plt.tight_layout(h_pad=0.75)
                if i % 50 == 0:
                    print(i, "of", length_of_movie,
                          "Making label_watershed_method video")
                # plt.tight_layout(h_pad=0.75)
                # plt.title(paramiters, y=1.08)
                writer.grab_frame()
                
        for c in ones_with_color:
            c[:,0,0] = 0
        


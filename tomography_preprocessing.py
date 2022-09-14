import tomopy
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy
from scipy.ndimage import median_filter
import dipy
from dipy.denoise.noise_estimate import estimate_sigma
from dipy.denoise.nlmeans import nlmeans
import os


class find_center:
    def __init__(self, fname):
        self.out = None
        self.arr = np.load(fname)
        print("Loaded ", fname)
        # the place the ring is usually centered if axis of rotation
        # is the same as the 0 axis of the volume array
        self.half_min_side_length = \
            min([self.arr.shape[1], self.arr.shape[2]]) / 2.0
        self.slice_index = self.arr.shape[0] // 2
        
        # gui variables
        self.locked = False
        self.finalized = False
        
        self.x = []
        self.y = []
        self.scale = []
    
    # this is usually good enough
    def use_center_of_image(self):
        self.out = np.array([int(self.half_min_side_length),
                             int(self.half_min_side_length)])
        return self.out
    
    # however if you see rings in the de-noised data
    # it might be worth selecting the center yourself
    def select_with_gui(self):
        
        # the size of the circles (this is arbitrary) just a visual aid.
        def get_s():
            return len(self.x) * self.half_min_side_length / 3.0 + 2
        
        def onclick(event):
            if self.finalized:
                return self.out
            
            s = get_s()
            
            if event.button == 1 and event.inaxes is not None:
                # guess a location (will place a circle )
                finalized = False
                print("Guess: x: %s y: %s" % (event.xdata, event.ydata))
                print('you have guessed this many times:', len(self.x) + 1)
                # left click adds a circle
                print(event)
                print(event.button)
                self.scale.append(np.pi * (s ** 2))
                self.x.append(event.xdata)
                self.y.append(event.ydata)
                # clear frame
                plt.clf()
                plt.imshow(self.arr[self.slice_index], cmap='gray',
                           origin='lower')
                plt.xlabel("e1")
                plt.ylabel("e2")
                plt.title("left click to add circles," +
                          " (recomend using at least 5)\n" +
                          "right click to remove last circle or approve mean")
                plt.scatter(self.x, self.y, self.scale, facecolors='none',
                            edgecolors='r')
                # inform matplotlib of the new data
                plt.draw()  # redraw
            
            if event.button == 3:
                # right clicking gives you the option to say your are done,
                # you want to remove a circle, or that it was a mistake
                done = {'d', 'done'}
                no = {'no', 'n'}
                remove = {'r', 'remove'}
                sys.stdout.write("\nEnter done and close window to lock in selected centers\n" +
                                 "remove to remove previouse circle\n" +
                                 "no to add more circles\n")
                choice = input().lower()
                if choice in remove:
                    print("Remove Last Click")
                    if len(self.x) >= 1:
                        self.x = self.x[:-1]
                        self.y = self.y[:-1]
                        self.scale = self.scale[:-1]
                    
                    # clear frame
                    plt.clf()
                    plt.imshow(self.arr[self.slice_index], cmap='gray',
                               origin='lower')
                    plt.xlabel("e1")
                    plt.ylabel("e2")
                    plt.title(
                        "left click to add circles, " +
                        "center will be average of the circles\n" +
                        " right click to remove last circle or approve mean")
                    plt.scatter(self.x, self.y, self.scale, facecolors='none',
                                edgecolors='r')
                    # inform matplotlib of the new data
                    plt.draw()  # redraw
                elif choice in no:
                    plt.clf()
                    plt.imshow(self.arr[self.slice_index], cmap='gray',
                               origin='lower')
                    plt.xlabel("e1")
                    plt.ylabel("e2")
                    plt.title("left click to add circles, " +
                              "(recomend using at least 5)" +
                              "\n right click to remove last" +
                              " circle or approve mean")
                    plt.scatter(self.x, self.y, self.scale, facecolors='none',
                                edgecolors='r')
                    # inform matplotlib of the new data
                    plt.draw()  # redraw
                elif choice in done:
                    self.finalized = True
                    
                    best_x = [np.mean(self.x), np.mean(self.x)]
                    best_y = [np.mean(self.y), np.mean(self.y)]
                    plt.clf()
                    plt.imshow(self.arr[self.slice_index],
                               cmap='gray', origin='lower')
                    plt.xlabel("e1")
                    plt.ylabel("e2")
                    plt.title("To approve, close the plot, else right click\
                     and select no")
                    plt.scatter(best_x, best_y, [20, 30], facecolors='none',
                                edgecolors='r')
                    # inform matplotlib of the new data
                    plt.draw()  # redraw
                
                else:
                    sys.stdout.write("Please respond with 'done'," +
                                     " 'remove', or 'no'\n")
        
        def motion_notify(event):
            # this places the purple ring and lets you place it where
            # it makes sense
            if self.locked or self.finalized:
                return
            s = get_s()  # scale for this particular instance
            s = np.pi * (s ** 2)
            # clear frame and plot new stuff this time
            # including a hypothetical circle that may or may not be added
            plt.clf()
            plt.imshow(self.arr[self.slice_index], cmap='gray', origin='lower')
            plt.xlabel("e1")
            plt.ylabel("e2")
            plt.title(
                "left click to add circles, \
                (recomend using at least 5)\n \
                 right click to remove last circle or approve mean")
            plt.scatter([event.xdata], [event.ydata], s, facecolors='none',
                        edgecolors='purple')  # adds candidate cirlce
            
            if len(self.scale) != 0:
                plt.scatter(self.x, self.y, self.scale, facecolors='none',
                            edgecolors='r')
                # inform matplotlib of the new data
            plt.draw()  # redraw
        
        def handle_close(event):
            if not self.finalized:
                sys.stdout.write("ERROR: Program not\
                 properly exited bad mean returned\n")
            
            elif len(self.x) == 0:
                sys.stdout.write("ERROR: No Circles have been selected,\
                 please select the center of the rings\n")
            
            else:
                sys.stdout.write("Good Exit, mean over " +
                                 "%s points e1: %s ,e2: %s\n" %
                                 (len(self.x), np.mean(self.x),
                                  np.mean(self.y)))
                self.locked = True
                self.out = np.array([np.mean(self.x), np.mean(self.y)])
                self.out = np.int32(np.around(self.out))
                return self.out
        
        fig, ax = plt.subplots()
        fig.canvas.mpl_connect('button_release_event', onclick)
        fig.canvas.mpl_connect('motion_notify_event', motion_notify)
        fig.canvas.mpl_connect('close_event', handle_close)
        
        plt.title("Line up circles with ring artifact to\n" +
                  "determine it's center for de noising" +
                  " left Click to begin")
        plt.xlabel("e1")
        plt.ylabel("e2")
        plt.show()
        # plt.draw()
        
        return self.out



class tomo_preprocessor:
    
    def __init__(self, outputfolder,
                 fname, center_x, center_y, save_steps_and_images=False,
                 mask_radius=None):
        self.debug = save_steps_and_images
        self.fname = outputfolder + fname
        # this is the volume to be denoised and deringed
        self.volume = np.load(fname)
        self.center_x = center_x
        self.center_y = center_y
        # if less than 2*width this will mask a cylinder
        # centered at the selected centerx center
        self.mask_radius = mask_radius

        #save starting volume to output folder
        if self.debug:
            np.save(fname, self.volume)

        # clip it so that the outer 10% is not gonna make stuff very wierd.
        self.vmin, self.vmax = scipy.stats.scoreatpercentile(self.volume, (0.5, 99.5))
        self.volume = np.clip(self.volume, self.vmin, self.vmax)
        smallest_val = np.amin(self.volume)
        self.volume = self.volume - smallest_val
        # this makes the smallest 0.5 values equal to zero

        if self.debug:
            np.save(self.fname[:-4] + "_center.npy", np.array([center_x, center_y]))
            # starting distribution display
            tmp = self.volume.copy()
            tmp = tmp.flatten()
            plt.hist(tmp, bins=100)
            title = self.fname[:-4]
            title = title.split("/")
            title = title[-1]
            plt.title(title +
                      "\nStarting Value distribution in Reconstruction")
            plt.savefig(self.fname[:-4] + "_starting_value_distribution.png")
            
        
        #mask if you have wierd stuff outside the packing sample
        self.mask_3d = None
        if self.mask_radius is not None:
            # make mask if selected
            mask_x = int(self.center_x)
            mask_y = int(self.center_y)
            dim1 = self.volume.shape[1]
            dim2 = self.volume.shape[2]
            y_grid, x_grid = np.ogrid[-mask_x:dim1 - mask_x,
                                      -mask_y:dim2 - mask_y]
            mask_2d = y_grid * y_grid + x_grid * x_grid <= self.mask_radius ** 2
            mask_2d = np.int8(mask_2d)
            self.mask_3d = np.int8(np.zeros(self.volume.shape))
            for i in np.arange(self.volume.shape[0]):
                self.mask_3d[i] = mask_2d[:]
        
    
    #one kind of filter that does work, other smoothing filters might also work
    #the purpose of this is to get rid of large differences
    def apply_non_local_means(self, pr=1, br=3):
        # first use non local means to de noise the volume
        
        dipy_sigma = estimate_sigma(self.volume, N=0)
        self.volume = dipy.denoise.nlmeans.nlmeans(self.volume[:],
                                                     mask=self.mask_3d,
                                                     sigma=dipy_sigma,
                                                     patch_radius=pr,
                                                     block_radius=br,
                                                     rician=False)
        
        #self.volume = median_filter(self.volume, size=3)
        return self.volume
    
    # de ring using tomopy
    def apply_de_ring(self, rwidth=9):
        # by trial and error 9 seems good as a defualt
        # (it has to do with how wide the rings can be)
        self.volume = tomopy.misc.corr.remove_ring(self.volume[:],
                                                     center_x=self.center_x,
                                                     center_y=self.center_y,
                                                     rwidth=rwidth)
        
        min_val = np.amin(self.volume)
        self.volume = self.volume - min_val
        
        if self.mask_radius is not None:
            self.volume[self.mask_3d == 0] = 0

        return self.volume
        
    def save_output(self):
        if self.debug:
            
            np.save(self.fname[:-4] + "_preprocessed.npy", self.volume)
    
            # for display final distribution
            tmp = self.volume.copy()
            tmp = tmp.flatten()
            plt.close()
            plt.hist(tmp, bins=100)
            plt.xlim(xmin=0.15)
            title = self.fname[:-4]
            title = title.split("/")
            title = title[-1]
            plt.title(title +
                      "\nFinal (denoised) Value distribution" +
                      " in Reconstruction, not including 0")
            plt.savefig(self.fname[:-4] + "_preprocessed_value_distribution.png")

        return self.volume

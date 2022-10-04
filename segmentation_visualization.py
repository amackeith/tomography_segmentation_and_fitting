import glob

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap
from watershed_segmentation import rand_cmap
from scipy import ndimage

### this is code to compare two different segmentations
def compare_segmentation(density_file, segments, nm, direc, rng=None):
    colors = rand_cmap(10000, verbose=False)
    segments[:, 0, 0] = 10000
    os.chdir(direc)

    seg_cmap = colors
    #seg_cmap[:, -1] = np.linspace(0, 1, seg_cmap.N)
    #seg_cmap = ListedColormap(seg_cmap)

    FFMpegWriter = animation.writers['ffmpeg']
    writer = FFMpegWriter(fps=5)

    fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    l, = plt.plot([], [], 'k-o')

    length_of_movie = segments.shape[2]

    print(segments.shape, density_file.shape)
    if rng == None:
        rng = [0, length_of_movie]
    #with writer.saving(fig, "hello" + ".mp4", length_of_movie):
    if True:
        for i in range(rng[0], rng[1]):

            ###### without this block of code this process will take like
            ###### 45 minutes
            ax.clear()
            fig, ax = plt.subplots(1, 1, figsize=(9, 9))

            ax.imshow(density_file[:,:,i],
                            cmap='gray', origin='lower')
            ax.set_title(f"Slice {i:03}")

            ax.imshow(segments[:,:,i],
                            origin='lower',
                            cmap=seg_cmap, interpolation='none', alpha=0.5)

            if i % 1 == 0:
                print(i, "of", length_of_movie,
                      "Making label_watershed_method video")
            # plt.tight_layout(h_pad=0.75)
            # plt.title(paramiters, y=1.08)

            #plt.show()
            plt.savefig(f"{nm}{i:03}.png")
            #writer.grab_frame()




def compare_single_cell(binary_desnity_file, density_file, segments, binary_single_cell, nm, direc, rng=None):
    os.chdir(direc)
    #labels_of_note = segments.copy()
    #single_cell  = binary_single_cell.copy()
    x = ndimage.center_of_mass(binary_single_cell)

    #labels_of_note[binary_single_cell == 0] = 0;
    #labels_of_note = list(set(labels_of_note.flatten()))
    all_labels = list(set(segments.flatten()))
    x = np.round(x)
    x = np.int32(x)
    print(x)
    labels_of_note = [segments[x[0], x[1], x[2]], ]
    visual_padding = 10
    # zero out all non notable labels
    print(labels_of_note)
    print(len(labels_of_note), len(all_labels))
    for i in all_labels:
        if i not in labels_of_note:
            segments[segments == i] = 0

    import scipy.io
    mat_dict = {'arr':segments[30:-30, 30:-30, 30:-30]}
    scipy.io.savemat("watershed_single_cell.mat", mat_dict)
    #exit()



    #labels_of_note = segments.copy()
    #labels_of_note[binary_single_cell == 0] = 0;
    #labels_of_note = list(set(labels_of_note.flatten()))
    #print(labels_of_note)

    #print(len(labels_of_note), len(all_labels))

    #xstack = np.any(binary_single_cell, axis=1)
    #ystack = np.any(binary_single_cell, axis=0)
    #zstack = np.any(binary_single_cell, axis=2)
    #ymin, ymax = np.where(ystack)[0][[0, -1]]
    #xmin, xmax = np.where(xstack)[0][[0, -1]]
    #zmin, zmax = np.where(zstack)[0][[0, -1]]


    #rng = [max(0, zmin - visual_padding), min(zmax + visual_padding, segments.size[2])]
    colors = rand_cmap(10000, verbose=False)

    colors2 = rand_cmap(10000, verbose=False, random_seed=88)
    segments[:, 0, 0] = 10000


    seg_cmap = colors
    #seg_cmap[:, -1] = np.linspace(0, 1, seg_cmap.N)
    #seg_cmap = ListedColormap(seg_cmap)

    FFMpegWriter = animation.writers['ffmpeg']
    writer = FFMpegWriter(fps=5)

    fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    l, = plt.plot([], [], 'k-o')

    length_of_movie = segments.shape[2]

    print(segments.shape, density_file.shape)
    if rng == None:
        rng = [0, length_of_movie]
    #with writer.saving(fig, "hello" + ".mp4", length_of_movie):
    if True:
        for i in range(rng[0], rng[1]):
            print("HI")
            ###### without this block of code this process will take like
            ###### 45 minutes
            plt.close('all')

            if True:
                fig, ax = plt.subplots(2, 3, figsize=(9, 9))
                for ll in range(2):
                    for mm in range(3):
                        ax[ll,mm].clear()

                ax[0, 0].set_title(f"Slice {i:03} - Density - handmade")
                ax[0, 0].imshow(density_file[:, :, i],
                                         cmap='gray')
                ax[0, 0].imshow(binary_single_cell[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.5)

                ax[0, 1].imshow(density_file[:, :, i],
                                cmap='gray')
                ax[0, 1].set_title("Density - watershed")
                ax[0, 1].imshow(segments[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.5)


                ax[0, 2].imshow(density_file[:, :, i],
                                cmap='gray')
                ax[0, 2].set_title("Density - both")
                ax[0, 2].imshow(segments[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.25)
                ax[0, 2].imshow(binary_single_cell[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.25)

                ax[1, 0].set_title(f"Slice {i:03} - Binary - handmade")
                ax[1, 0].imshow(binary_desnity_file[:, :, i],
                                         cmap='gray')
                ax[1, 0].imshow(binary_single_cell[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.5)

                ax[1, 1].imshow(binary_desnity_file[:, :, i],
                                cmap='gray')
                ax[1, 1].set_title("Binary - watershed")
                ax[1, 1].imshow(segments[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.5)


                ax[1, 2].imshow(binary_desnity_file[:, :, i],
                                cmap='gray')
                ax[1, 2].set_title("Binary - both")
                ax[1, 2].imshow(segments[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.25)
                ax[1, 2].imshow(binary_single_cell[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.25)

                plt.savefig(f"{nm}_combined_{i:03}.png")
                # writer.grab_frame()

            if True:

                fig, ax = plt.subplots(1, 1, figsize=(9, 9))

                ax.set_title(f"Slice {i:03} - Density - handmade")
                ax.imshow(density_file[:, :, i],
                                cmap='gray')
                ax.imshow(binary_single_cell[:, :, i],
                                cmap=colors2, interpolation='none', alpha=0.5)
                plt.tight_layout()
                plt.savefig(f"{nm}_density_handmade_{i:03}.png")
                fig, ax = plt.subplots(1, 1, figsize=(9, 9))

                ax.imshow(density_file[:, :, i],
                                cmap='gray')
                ax.set_title(f"Slice {i:03} Density - watershed")
                ax.imshow(segments[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.5)

                plt.tight_layout()
                plt.savefig(f"{nm}_density_watershed_{i:03}.png")
                fig, ax = plt.subplots(1, 1, figsize=(9, 9))

                ax.set_title(f"Slice {i:03} - Binary - handmade")
                ax.imshow(binary_desnity_file[:, :, i],
                                cmap='gray')
                ax.imshow(binary_single_cell[:, :, i],
                                cmap=colors2, interpolation='none', alpha=0.5)

                plt.tight_layout()
                plt.savefig(f"{nm}_binary_handmade_{i:03}.png")
                fig, ax = plt.subplots(1, 1, figsize=(9, 9))

                ax.imshow(binary_desnity_file[:, :, i],
                                cmap='gray')
                ax.set_title(f"Slice {i:03} Binary - watershed")
                ax.imshow(segments[:, :, i],
                                cmap=seg_cmap, interpolation='none', alpha=0.5)

                plt.tight_layout()
                plt.savefig(f"{nm}_binary_watershed_{i:03}.png")

                fig, ax = plt.subplots(1, 1, figsize=(9, 9))

                ax.set_title(f"Slice {i:03} Density - both")
                ax.imshow(density_file[:, :, i],
                          cmap='gray')
                ax.imshow(segments[:, :, i],
                          cmap=seg_cmap, interpolation='none', alpha=0.25)
                ax.imshow(binary_single_cell[:, :, i],
                          cmap=colors2, interpolation='none', alpha=0.25)

                plt.tight_layout()
                plt.savefig(f"{nm}_density_both_{i:03}.png")

                fig, ax = plt.subplots(1, 1, figsize=(9, 9))

                ax.imshow(binary_desnity_file[:, :, i],
                                cmap='gray')
                ax.set_title(f"Slice {i:03} Binary - both")
                ax.imshow(segments[:, :, i],
                          cmap=seg_cmap, interpolation='none', alpha=0.25)
                ax.imshow(binary_single_cell[:, :, i],
                          cmap=colors2, interpolation='none', alpha=0.25)
                plt.tight_layout()
                plt.savefig(f"{nm}_binary_both_{i:03}.png")






            if i % 1 == 0:
                print(i, "of", rng[1],
                      "Making label_watershed_method video")
            # plt.tight_layout(h_pad=0.75)
            # plt.title(paramiters, y=1.08)

            #plt.show()
            #writer.grab_frame()

        names = ["_binary_both_", "_density_both_", "_binary_handmade_", "_binary_watershed_", "_density_watershed_", "_density_handmade_",
                 "_combined_"]

        for name in names:
            bashCommand = f"convert {direc}{nm}{name}*.png {direc}{nm}{name[:-1]}.mp4"
            os.system(bashCommand)



lst = glob.glob("~/Documents/ohern/experimental_leaf_data/jscheelChunks_20220919/*padded.npy")

human_made_sm = np.load("/home/arthur/Documents/ohern/experimental_leaf_data/jscheelChunks_20220919/cell_0157_padded.npy")
binary_chunks = np.load("/home/arthur/Documents/ohern/experimental_leaf_data/jscheelChunks_20220919/chunk_cell_0157_padded.npy")
segmented = np.load("/home/arthur/Documents/ohern/experimental_leaf_data/jscheelChunks_20220919/Chunk_0157_99.0_0.0/chunk_cell_0157_padded0.5_second_grow_labels.npy")
grey_file = np.load("/home/arthur/Documents/ohern/experimental_leaf_data/jscheelChunks_20220919/chunk_cell_0157_RGB_padded.npy")

sm_shape = human_made_sm.shape
pos = [109, 143, 160]
pos = [97,   183,   264]
pos = [90, 145, 165];

padding = 30
pos = np.array(pos)
human_made_sm = np.rot90(human_made_sm, 2, (0,2));
human_made = np.zeros(binary_chunks.shape)
human_made[pos[0]:pos[0]+sm_shape[0], pos[1]:pos[1]+sm_shape[1], pos[2]:pos[2]+sm_shape[2]] = human_made_sm
#109   143   289
print(binary_chunks.shape, human_made.shape)


#binary_chunks = np.swapaxes(binary_chunks, 2,0);
#human_made = np.swapaxes(human_made, 2,0);
#grey_file = np.swapaxes(grey_file, 2,0);
#segmented = np.swapaxes(segmented, 2,0);



compare_single_cell(binary_chunks, grey_file, segmented, human_made, 'chunk_0157_00',"/home/arthur/Documents/ohern/experimental_leaf_data/jscheelChunks_20220919/Chunk_0157_99.0_0.0/", rng=[pos[2], pos[2]+human_made_sm.shape[2]])
exit()


if __name__ == "__main__":
    #direc = "/home/arthur/Documents/ohern/experimental_leaf_data/analysis_code/chunk_0096_padded2/"
    direc = os.path.abspath(sys.argv[1]) + "/"
    # python /home/arthur/Documents/coding_projects/tomography_segmentation_and_fitting/segmentation_visualization.py Chunk_0157_99.75_100.0/
    nm = direc.split("/")
    nm = nm[-2]
    print(nm)
    density_fname = glob.glob(direc + "*original.npy")[0]
    labels_fname = glob.glob(direc + "*second_grow_labels.npy")[0]
    density = np.load(density_fname, mmap_mode='r')
    labels = np.load(labels_fname, mmap_mode='r')
    labels = labels * 1
    density = density
    compare_segmentation(density, labels, direc, nm)

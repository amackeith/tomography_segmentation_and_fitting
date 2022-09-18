import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap
from watershed_segmentation import rand_cmap


def compare_segmentation(density_file, segments):
    colors = rand_cmap(10000, verbose=False)
    segments[:, 0, 0] = 10000


    seg_cmap = colors
    #seg_cmap[:, -1] = np.linspace(0, 1, seg_cmap.N)
    #seg_cmap = ListedColormap(seg_cmap)

    FFMpegWriter = animation.writers['ffmpeg']
    writer = FFMpegWriter(fps=5)

    fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    l, = plt.plot([], [], 'k-o')

    length_of_movie = segments.shape[0]

    print(segments.shape, density_file.shape)

    #with writer.saving(fig, "hello" + ".mp4", length_of_movie):
    if True:
        for i in range(0, length_of_movie):

            ###### without this block of code this process will take like
            ###### 45 minutes
            ax.clear()
            fig, ax = plt.subplots(1, 1, figsize=(9, 9))

            ax.imshow(density_file[i],
                            cmap='gray', origin='lower')
            ax.set_title(f"Slice {i}")

            ax.imshow(segments[i],
                            origin='lower',
                            cmap=seg_cmap, interpolation='none', alpha=0.5)

            if i % 1 == 0:
                print(i, "of", length_of_movie,
                      "Making label_watershed_method video")
            # plt.tight_layout(h_pad=0.75)
            # plt.title(paramiters, y=1.08)

            #plt.show()
            plt.savefig(f"b{i}.png")
            #writer.grab_frame()




if __name__ == "__main__":
    direc = "/home/arthur/Documents/ohern/experimental_leaf_data/analysis_code/chunk_0096_padded2/"
    density = np.load(direc + "chunk_0096_padded0.5_binary.npy", mmap_mode='r')
    labels = np.load(direc + "chunk_0096_padded0.5_first_grow_labels.npy", mmap_mode='r')
    labels = labels * 1
    density = density
    compare_segmentation(density, labels)
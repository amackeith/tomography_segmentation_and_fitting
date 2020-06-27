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


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
import sys

if len(sys.argv) != 2:
	print("binary_threshold_selection.py not passed the right number of args")
	print("Should be binary_threshold_selection.py filename")
	exit()
else:
	fname = sys.argv[1]
	np_image = np.load(fname)



histvals = np_image.copy()
histvals = histvals.flatten()


fig, ax = plt.subplots(1,3, figsize=(16,8))
ax[0].hist(histvals, bins=100)
(n,bins,patches) = plt.hist(histvals, bins=100)
ymax = np.max(n[10:])*1.05
ax[0].set_xlim(xmin=0.05)
ax[0].set_ylim(ymax=ymax)
ax[0].set_title("Tomography value histogram (not including 0)")
fig.subplots_adjust(left=0.1, bottom=0.25)
min0 = 0
max0 = (np_image.shape[0]-1)/2
im1  = ax[1].imshow(np_image[max0])
ax[1].set_title("Select Threshold using slider")





example = np.load("example_binary_thresh.npy")
min1 = 0
max1 = (example.shape[0]-1)/2
im2  = ax[2].imshow(example[max1])




axcolor = 'lightgoldenrodyellow'
axmin = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axmax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_example = fig.add_axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)

half = (np.min(bins)+ np.max(bins))/2
threshold_slider = Slider(axmin, 'Threshold', np.min(bins), np.max(bins), valinit=half)
slice_slider = Slider(axmax, 'Slice', 0, np_image.shape[0]-1, valinit=max0)
example_slider = Slider(ax_example, 'Example Slice', 0, example.shape[0]-1, valinit=max1)


def update(val):
	ax[1].cla()
	im1 = ax[1].imshow(np.array((np_image>threshold_slider.val)*1.0)[int(slice_slider.val)])
	ax[2].cla()
	im2 = ax[2].imshow(example[int(example_slider.val)])
	fig.canvas.draw()

threshold_slider.on_changed(update)
slice_slider.on_changed(update)
example_slider.on_changed(update)
plt.show()

thresh = np.array([threshold_slider.val])
print(thresh.shape)
print(thresh)
np.save("%s_selected_binary_threshold.npy" % fname, thresh)

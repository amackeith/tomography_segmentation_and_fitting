
import tomopy
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy
import dipy
from dipy.denoise.noise_estimate import estimate_sigma
from dipy.denoise.nlmeans import nlmeans
import os




if len(sys.argv) != 3:
	print("make_de_noised_de_ringed_given_cenxy.py passed wrong number of args")
	print("should be: python make_de_noised_de_ringed_given_cenxy volume_file_name.npy")
	exit()
elif not os.path.exists("center_of_ring.npy"):
	print("center_of_ring.npy doesn't exists, check for issues with select center")
	exit()
else:
	fname = sys.argv[1]
	max_r = int(float(sys.argv[2]))
	cen = np.load("center_of_ring.npy")
	center_x = cen[0]
	center_y = cen[1]

#this is the volume to be denoised and deringed
test_vol = np.load(fname)


#clip it so that the outer 10% is not gonna make stuff very wierd.
vmin, vmax = scipy.stats.scoreatpercentile(test_vol, (0.5, 99.5))
test_vol = np.clip(test_vol, vmin, vmax)
smallest_val = np.amin(test_vol)
test_vol = test_vol - smallest_val #this makes the smallest 0.5 values equal to zero

#starting distribution display
tmp = test_vol.copy()
tmp = tmp.flatten()
plt.hist(tmp, bins=100)
plt.title("Starting Value distribution in Reconstruction")
plt.savefig("starting_value_distribution.png")


#make mask
mask_x = int(center_x)
mask_y = int(center_y)
dim1 = test_vol.shape[1]
dim2 = test_vol.shape[2]
y_grid, x_grid = np.ogrid[-mask_x:dim1 - mask_x, -mask_y: dim2 - mask_y]
mask_2d = y_grid*y_grid + x_grid*x_grid <= max_r**2
mask_2d = np.int8(mask_2d)
mask_3d = np.int8(np.zeros(test_vol.shape))
for i in np.arange(test_vol.shape[0]):
	mask_3d[i] = mask_2d[:]



#first use non local means to de noise the volume
dipy_sigma = estimate_sigma(test_vol, N=0)
test_vol = dipy.denoise.nlmeans.nlmeans(test_vol[:], mask=mask_3d,sigma=dipy_sigma, patch_radius= 1, block_radius = 3, rician=False)




#de ring using tomopy
rwidth = 9
#9 is lokking p good (it has to do with how wide the rings can be)
test_vol = tomopy.misc.corr.remove_ring(test_vol[:], center_x=center_x, center_y=center_y, rwidth=rwidth)


#shift everything so that it is positive, then re-mask so that the background is still 0
min_val = np.amin(test_vol)
test_vol = test_vol-min_val
test_vol[mask_3d==0]=0
#save file
np.save("%s_de_noised_de_ringed.npy" % fname, test_vol)


#for display final distribution 
tmp = test_vol.copy()
tmp = tmp.flatten()
plt.close()
plt.hist(tmp, bins=100)
plt.xlim(xmin=0.15)
plt.title("Final (denoised) Value distribution in Reconstruction, not including 0")
plt.savefig("final_denoisded_value_distribution.png")




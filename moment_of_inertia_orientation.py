'''
Author: Arthur MacKeith
Date: August 20, 2018

Based on methods of Fabian M. Schaller, Matthias Schoter et al.
"Tomographic Analysis of Jammed Ellipsiod Packings"
in "Powders and Grains 2013"

This code finds the orientation of uniaxial particles by finding thier largest moment of inertia eigen value
and saying that the vector coresponding is the principle axis. A slight edit will have to be made for "prolate" shapes
see tag #edit_here in moment_of_inertia_processing_multi_core.py
'''
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import sys


cpp_enabled = True
try:
    import watershedtools_cpp.watershedtools_cpp as watershedtools_cpp
except:
    print("CPP code not enabled")
    print("to build and enable it run ./build in the watershedtools_cpp folder")
    print("Proceeding using SLOWER python code")
    cpp_enabled = False


#in arr is a binary lentil

def find_principal_axis(in_arr):
    arr = in_arr
    #com will be used as origin
    #print "COM:",com
    
    

    moment_of_inertia = np.zeros((3,3))*1.0

    if cpp_enabled:
        com = np.zeros(3) * 1.0
        watershedtools_cpp.center_of_mass(arr, com)
        
        watershedtools_cpp.calculate_moment_of_inertia(in_arr,
                                                       moment_of_inertia,
                                                       com[0], com[1], com[2])
    else:
        com = np.array(ndi.measurements.center_of_mass(arr))
        shp = in_arr.shape
        for i in range(shp[0]):
            for j in range(shp[1]):
                for k in range(shp[2]):
                    if arr[i,j,k]==1:
                        #x,y,z is the distance for the origin in each direction
                        x,y,z = np.array([i,j,k])-com
                        moment_of_inertia[0,0] =\
                            moment_of_inertia[0,0]+y**2+z**2	#I_xx
                        moment_of_inertia[1,1] =\
                            moment_of_inertia[1,1]+x**2+z**2	#I_yy
                        moment_of_inertia[2,2] =\
                            moment_of_inertia[2,2]+x**2+y**2 #I_zz
                        moment_of_inertia[0,1] =\
                            moment_of_inertia[0,1]-x*y #I_xy
                        moment_of_inertia[1,2] =\
                            moment_of_inertia[1,2]-y*z	#I_yz
                        moment_of_inertia[0,2] =\
                            moment_of_inertia[0,2]-x*z	#I_xz
    
        moment_of_inertia[1,0] = moment_of_inertia[0,1] #I_yx
        moment_of_inertia[2,1] = moment_of_inertia[1,2]	#I_zy
        moment_of_inertia[2,0] = moment_of_inertia[0,2]	#I_zx
        #
    #such that eigen_values[i] is the ith Eval,
    #eigen_vectors[:,i] is the coresponding eigen vector
    eigen_values, eigen_vectors = np.linalg.eig(moment_of_inertia)


    #print "EVAL", eigen_values
    #print "EVEC\n",eigen_vectors


    #ident = np.array([[1.0,0,0], [0,1.0,0], [0,0,1.0]])

    #proof
    #print "should be zero", np.linalg.det(moment_of_inertia-ident*eigen_values[0]),np.linalg.det(moment_of_inertia-ident*eigen_values[1]),np.linalg.det(moment_of_inertia-ident*eigen_values[2])
    #print moment_of_inertia
    return eigen_values, eigen_vectors, com


#this is sus aff rn
def get_phi_and_theta_from_principle_axis(pi_vec):
    #since theta in 0,pi and phi in 0, pi caputres every orientation
    #and since in those ranges r has a positive x component
    #if x dot r is < 0 use -r instead
    x_unit = np.array([1,0,0])
    z_unit = np.array([0,0,1])
    neg_y = np.array([0,-1,0])

    if np.dot(x_unit,pi_vec)<0:
        pi_vec = -1.0*pi_vec

    #these projections use the euler angles to figure out what phi and theta are

    r_unit_vec = pi_vec/np.linalg.norm(pi_vec)
    theta = np.arccos(np.dot(z_unit, r_unit_vec))
    r_proj_on_xy = r_unit_vec.copy()
    r_proj_on_xy[2] = 0
    if np.linalg.norm(r_proj_on_xy)==0.0:
        phi=0.0
    else:
        r_proj_on_xy = r_proj_on_xy/np.linalg.norm(r_proj_on_xy)
        phi = np.arccos(np.dot(neg_y, r_proj_on_xy))
    
    return (phi, theta)





'''
The rest of the code is for testing and may be a useful guide if you need to benchmark this for accuracy





padding = 30 ## voxels on each side of xy
xvals, yvals, zvals = np.meshgrid(range(2*padding),range(2*padding),range(2*padding),indexing="ij")
radius = 22.8 ## size scale of lentils
def make_rot_mat_simplified(phi, theta):
    rotmat = np.zeros((3,3))
    rotmat[0,0] = np.cos(phi)
    rotmat[0,1] = np.sin(phi)
    rotmat[0,2] = 0.

    rotmat[1,0] = -np.cos(theta)*np.sin(phi)
    rotmat[1,1] = np.cos(theta)*np.cos(phi)
    rotmat[1,2] = np.sin(theta)

    rotmat[2,0] = np.sin(theta)*np.sin(phi)
    rotmat[2,1] = -np.sin(theta)*np.cos(phi)
    rotmat[2,2] = np.cos(theta)

    return rotmat

def create_lentil_phi_theta(params):
    phi,theta = params
    phi*=np.pi; theta*=np.pi
    ## in_arr is a bunch of zeros, xyz is deviations in voxels from the center of in_arr
    in_arr = np.zeros((2*padding,2*padding,2*padding))
    cenx = padding; ceny = padding; cenz = padding
    ## a lentil is defined by two points -- start off with the radii of these spheres being r and their separation is set by the angle of the lentil
    ## 45 deg lentils have the two points separated by sqrt2 * r
    ## the two points are initially aligned with the z axis
    ## do the rotation first, so that the points start out at (0,0,+- r sqrt2 / 2), then rotate according to
    rot_mat = make_rot_mat_simplified(phi, theta).T
    sphere1 = np.float32([0,0,radius*np.sqrt(2)/2.])
    sphere2 = np.float32([0,0,-radius*np.sqrt(2)/2.])

    sphere1 = np.dot(rot_mat,sphere1)
    sphere2 = np.dot(rot_mat,sphere2)
    on_vox = (((xvals-cenx-sphere1[0])**2+(yvals-ceny-sphere1[1])**2+(zvals-cenz-sphere1[2])**2)<radius**2) & \
                (((xvals-cenx-sphere2[0])**2+(yvals-ceny-sphere2[1])**2+(zvals-cenz-sphere2[2])**2)<radius**2)
    in_arr[on_vox] = 1
    return np.int8(in_arr)

def create_lentil(params):
    x,y,z,phi,theta = params
    phi*=np.pi; theta*=np.pi
    ## in_arr is a bunch of zeros, xyz is deviations in voxels from the center of in_arr
    in_arr = np.zeros((2*padding,2*padding,2*padding))
    cenx = padding; ceny = padding; cenz = padding
    ## a lentil is defined by two points -- start off with the radii of these spheres being r and their separation is set by the angle of the lentil
    ## 45 deg lentils have the two points separated by sqrt2 * r
    ## the two points are initially aligned with the z axis
    ## do the rotation first, so that the points start out at (0,0,+- r sqrt2 / 2), then rotate according to
    rot_mat = make_rot_mat_simplified(phi, theta).T
    sphere1 = np.float32([0,0,radius*np.sqrt(2)/2.])
    sphere2 = np.float32([0,0,-radius*np.sqrt(2)/2.])

    sphere1 = np.dot(rot_mat,sphere1)
    sphere2 = np.dot(rot_mat,sphere2)
    on_vox = (((xvals-cenx-x-sphere1[0])**2+(yvals-ceny-y-sphere1[1])**2+(zvals-cenz-z-sphere1[2])**2)<radius**2) & \
                (((xvals-cenx-x-sphere2[0])**2+(yvals-ceny-y-sphere2[1])**2+(zvals-cenz-z-sphere2[2])**2)<radius**2)
    in_arr[on_vox] = 1
    return np.int8(in_arr)


#this creates a lentil and places it in the space array
def place_and_create_lentil(params):
    x,y,z,phi,theta = params
    
    ## in_arr is a bunch of zeros, xyz is deviations in voxels from the center of in_arr
    in_arr = np.zeros((2*padding,2*padding,2*padding))
    cenx = padding; ceny = padding; cenz = padding
    ## a lentil is defined by two points -- start off with the radii of these spheres being r and their separation is set by the angle of the lentil
    ## 45 deg lentils have the two points separated by sqrt2 * r
    ## the two points are initially aligned with the z axis
    ## do the rotation first, so that the points start out at (0,0,+- r sqrt2 / 2), then rotate according to
    rot_mat = make_rot_mat_simplified(phi, theta).T
    sphere1 = np.float32([0,0,radius*np.sqrt(2)/2.])
    sphere2 = np.float32([0,0,-radius*np.sqrt(2)/2.])

    sphere1 = np.dot(rot_mat,sphere1)
    sphere2 = np.dot(rot_mat,sphere2)

    on_vox = (((xvals-cenx-x-sphere1[0])**2+(yvals-ceny-y-sphere1[1])**2+(zvals-cenz-z-sphere1[2])**2)<radius**2) & \
                (((xvals-cenx-x-sphere2[0])**2+(yvals-ceny-y-sphere2[1])**2+(zvals-cenz-z-sphere2[2])**2)<radius**2)

    in_arr[on_vox] = 1
    return np.int8(in_arr)






data = np.load(#This should be a file with labeled segments in it)
shp = np.array(data.shape)+100
nest = np.zeros(shp)
nest[50:50+data.shape[0],50:50+data.shape[1],50:50+data.shape[2]] = data
fitting = nest






lbls = fitting.copy()
lbls = lbls.flatten()
lbls = list(set(lbls))
lbls.sort()
vol = []

for particle_index in lbls[1:]:
    mask = np.array((fitting==particle_index)*1)

    mask_vol = np.sum(mask)
    if mask_vol>4:
        vol.append(mask_vol)
        


    np_where_info = np.where(mask==1)
    x_max = max(np_where_info[0])
    x_min = min(np_where_info[0])
    y_max = max(np_where_info[1])
    y_min = min(np_where_info[1])
    z_max = max(np_where_info[2])
    z_min = min(np_where_info[2])

    x_avg = (x_min+x_max)/2
    y_avg = (y_min+y_max)/2
    z_avg = (z_min+z_max)/2

    #mother = np.array(mask)
    #mother = np.int8((mother==1)*1) #this means there will only be on voxels when the two intersect
    
    arr = fitting[x_avg-padding:x_avg+padding, y_avg-padding:y_avg+padding, z_avg-padding:z_avg+padding]

    val,vec = find_principal_axis(arr)
    index_of_max_eigen_val = np.argmax(val)
    principle_axis = vec[:,index_of_max_eigen_val]




    #lent = mask[x_avg-padding:x_avg+padding, y_avg-padding:y_avg+padding, z_avg-padding:z_avg+padding]
    
    #arr = np.complex(arr)

    if mask_vol<1800 or (not np.array_equal(arr.shape, np.array([50,50,50]))):
        continue

    #mlab.contour3d(arr)

    print mask_vol, "label ", particle_index, "xyz ", x_avg-50, y_avg-50, z_avg-50
    
    mlab.pipeline.volume(mlab.pipeline.scalar_field(arr), color=(1,0,0))
    mlab.pipeline.volume(mlab.pipeline.scalar_field(lent), color=(0,1,0))
    mlab.show()

    # combine the objects into a single boolean array






x,y,z = np.indices((50,50,50))
cube1 = (x<4)&(y<4)&(z<4)

errors = []
for i in range(200):
    phi = np.random.random_sample()
    theta = np.random.random_sample()
    x0 = 5*np.random.random_sample()
    y0 = 5*np.random.random_sample()
    z0 = 5*np.random.random_sample()
    
    voxels = create_lentil((x0,y0,z0,phi,theta))
    com = ndi.measurements.center_of_mass(voxels)

    dff = np.array([x0+padding,y0+padding,z0+padding])-com
    print "COm diff", np.linalg.norm(dff)
    val,vec = find_principal_axis(voxels)
    index_of_max_eigen_val = np.argmax(val)
    principle_axis = vec[:,index_of_max_eigen_val]



    p, t = get_phi_and_theta_from_principle_axis(principle_axis)
    


    zhat = np.array([0,0,1])
    phi, theta, p, t = (phi*np.pi,theta*np.pi,p,t)
    found_pi = np.dot(make_rot_mat_simplified(p,t), zhat)
    actual_pi = np.dot(make_rot_mat_simplified(phi,theta), zhat)



    
    #print np.around([phi-p,theta-t], decimals=5)
    #print "error in com", dff
    #print "Angle (degrees) between actual pi and found pi:", np.rad2deg(np.arccos(np.dot(found_pi,actual_pi)))
    if np.rad2deg(np.arccos(np.dot(found_pi,actual_pi)))>0.2:
        print np.around([phi-p,theta-t], decimals=5), phi, p, theta, t
        print "error in com", dff
        print "Angle (degrees) between actual pi and found pi:", np.rad2deg(np.arccos(np.dot(found_pi,actual_pi)))

    errors.append(np.rad2deg(np.arccos(np.dot(found_pi,actual_pi))))

plt.hist(errors)
plt.xlabel("Angle between input and output principle axis, Degrees")
plt.savefig("error_in_orientation")

'''





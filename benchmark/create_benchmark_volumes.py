import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Queue
import multiprocessing as mp
import numpy as np
import scipy.ndimage as ndi

radius = 25.0
angle = 45.0
# angle of lentils
eta = np.deg2rad(angle)
padding = 200  ## voxels on each side of xy
# making lentil array functions
spacer = int(radius * 1.75)

xvals, yvals, zvals = np.meshgrid(range(2 * padding), range(2 * padding), range(2 * padding), indexing="ij")


def make_rot_mat_simplified(phi, theta):
    rotmat = np.zeros((3, 3))
    rotmat[0, 0] = np.cos(phi)
    rotmat[0, 1] = np.sin(phi)
    rotmat[0, 2] = 0.
    
    rotmat[1, 0] = -np.cos(theta) * np.sin(phi)
    rotmat[1, 1] = np.cos(theta) * np.cos(phi)
    rotmat[1, 2] = np.sin(theta)
    
    rotmat[2, 0] = np.sin(theta) * np.sin(phi)
    rotmat[2, 1] = -np.sin(theta) * np.cos(phi)
    rotmat[2, 2] = np.cos(theta)
    
    return rotmat


def place_lentil(in_arr, params):
    x, y, z, phi, theta = params
    phi *= np.pi;
    theta *= np.pi;
    # in_arr is a bunch of zeros, xyz is deviations in voxels from the center of in_arr
    cenx = padding;
    ceny = padding;
    cenz = padding
    # a lentil is defined by two points -- start off
    # with the radii of these spheres being r and their separation is set by the angle of the lentil
    # 45 deg lentils have the two points separated by sqrt2 * r
    # the two points are initially aligned with the z axis
    # do the rotation first, so that the points start out at (0,0,+- r sqrt2 / 2), then rotate according to
    rot_mat = make_rot_mat_simplified(phi, theta).T
    
    r = (radius) * np.cos(eta)  # np.sqrt(2.0)/2.0
    
    norm = np.float32([0, 0, 1])
    sphere1 = np.float32([0, 0, r])
    sphere2 = np.float32([0, 0, -r])
    
    norm = np.dot(rot_mat, norm)
    sphere1 = np.dot(rot_mat, sphere1)
    sphere2 = np.dot(rot_mat, sphere2)
    
    on_vox = (((xvals - cenx - x - sphere1[0]) ** 2 + (yvals - ceny - y - sphere1[1]) ** 2 + (
            zvals - cenz - z - sphere1[2]) ** 2) < radius ** 2) & \
             (((xvals - cenx - x - sphere2[0]) ** 2 + (yvals - ceny - y - sphere2[1]) ** 2 + (
                     zvals - cenz - z - sphere2[2]) ** 2) < radius ** 2)
    
    in_arr[on_vox] = 1
    
    return in_arr, norm


def grid_lattice():
    # lenses placed in grid (all w/ same orientation).
    np.random.seed(100)
    list_of_loc_angle = []
    output_array = np.zeros((2 * padding, 2 * padding, 2 * padding))
    
    for i in range(-3,4):
        print('grid lattice', i)
        for j in range(-3, 4):
            for k in range(-3, 4):
                # print(i, j, k)
                _, norm = place_lentil(output_array,
                                [i * spacer, j * spacer, k * spacer, 0, 0])
                list_of_loc_angle.append(
                    np.array([
                        np.array([i * spacer, j * spacer, k * spacer])+padding,
                        np.array([0, 0, 0]), norm]))
                

    return output_array, np.array(list_of_loc_angle)




def grid_lattice_random_orientations():
    # orientations will be random.
    np.random.seed(100)
    list_of_loc_angle = []
    output_array = np.zeros((2 * padding, 2 * padding, 2 * padding))
    for i in range(-3, 4):
        print('grid lat random orientation', i)
        for j in range(-3, 4):
            for k in range(-3, 4):
                # print(i, j, k)
                a = np.random.random()
                b = np.random.random()
                phi = 2 * np.pi * a
                theta = np.arccos(2 * b - 1)
                
                
                _, norm = place_lentil(output_array,
                                       [i * spacer, j * spacer, k * spacer, phi, theta])
                
                entry = np.array([np.array([i * spacer, j * spacer, k * spacer]) + padding,
                                  np.array([phi, theta, 0]), norm])
                list_of_loc_angle.append(
                    entry)
    
    return output_array, np.array(list_of_loc_angle)


def grid_random():
    # orientations will be random.
    np.random.seed(100)
    list_of_loc_angle = []
    output_array = np.zeros((2 * padding, 2 * padding, 2 * padding))
    
    for i in range(-3, 4):
        print('grid random orientation and position', i)
        for j in range(-3, 4):
            for k in range(-3, 4):
                a = np.random.random()
                b = np.random.random()
                phi = 2 * np.pi * a
                theta = np.arccos(2 * b - 1)
                noise = np.random.random(3) * 2.0
                pos = np.array([i * spacer, j * spacer, k * spacer]) + noise
                _, norm = place_lentil(output_array,
                                       [pos[0], pos[1], pos[2], phi, theta])
                
                list_of_loc_angle.append(
                    np.array([pos + padding,
                              np.array([phi, theta, 0]), norm]))
    
    return output_array, np.array(list_of_loc_angle)


# makes a lattice (with aligned orientations)
def make_lattice(que):
    a, a_la = grid_lattice()
    np.save('grid_lattice.npy', a)
    np.save('grid_lattice_la.npy', a_la)
    
    que.put(None)


# positions on grid, orientations randomized
def make_lattice_random_orientations(que):
    b, b_la = grid_lattice_random_orientations()
    np.save('grid_lattice_random_orientations.npy', b)
    np.save('grid_lattice_random_orientations_la.npy', b_la)
    print(np.sum(np.int8(b == 2)))
    
    que.put(None)


# random positions and orientations (though no touching)
def make_random(que):
    c, c_la = grid_random()
    np.save('grid_random.npy', c)
    np.save('grid_random_la.npy', c_la)
    print(np.sum(np.int8(c == 2)))
    
    que.put(None)


def main():
    q = Queue()
    x = Process(target=make_lattice, args=(q,))
    x.start()
    y = Process(target=make_lattice_random_orientations, args=(q,))
    y.start()
    z = Process(target=make_random, args=(q,))
    z.start()
    
    
    
    x.join()
    y.join()
    z.join()


if __name__ == '__main__':
    main()

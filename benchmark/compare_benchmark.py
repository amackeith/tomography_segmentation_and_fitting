import numpy as np
import matplotlib.pyplot as plt


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
    # rot_mat = make_rot_mat(phi, theta, psi).T
    # norms = (np.dot(rot_mat, np.float32([1,0,0])), np.dot(rot_mat, np.float32([0,1,0])), np.dot(rot_mat, np.float32([0,0,1])))
    norm1 = np.dot(rot_mat, np.float32([0, 0, 1]))  # because you can specify  the only normal in


def test_grid():
    true_values = np.load('grid_lattice_la.npy', allow_pickle=True)
    fit_values = np.load('grid_lattice/grid_lattice0.5_positions_orientation.npy',
                         allow_pickle=True)
    
    angle_diff = []
    for i in true_values:
        distance = 1000
        match = None
        for j in fit_values:
            x = np.sum(np.abs(i[0] - j[0]))
            
            if x < distance:
                distance = x
                match = j
        
        if distance != 0:
            print("Failed to match every test particle", distance)
            print(i)
            print(match)
            
            
        else:
            angle_diff.append(
                np.rad2deg(np.arccos(np.abs(np.dot(i[2], match[2])))))
    
    if len(angle_diff) != len(true_values):
        print("grid lattice FAILED, did not match this many particles: ",
              len(true_values) - len(angle_diff))
    else:
        print("Matched Position of Every Test Particle on lattice")


def test_random_orientations():
    true_values = np.load('grid_lattice_random_orientations_la.npy', allow_pickle=True)
    fit_values = np.load(
        'grid_lattice_random_orientations/grid_lattice_random_orientations0.5_positions_orientation.npy',
        allow_pickle=True)
    
    angle_diff = []
    for i in true_values:
        distance = 1000
        match = None
        for j in fit_values:
            x = np.sum(np.abs(i[0] - j[0]))
            
            if x < distance:
                distance = x
                match = j
        
        if distance != 0:
            print("Failed to match every test particle", distance)
            
        else:
            angle_diff.append(
                np.rad2deg(np.arccos(np.abs(np.dot(i[2], match[2])))))
            

    if len(angle_diff) != len(true_values):
        print("random orientation FAILED, did not match this many particles: ",
              len(true_values) - len(angle_diff))
        exit()
    
    
    plt.hist(angle_diff, bins=30)
    plt.xlabel("Degrees")
    plt.ylabel("Count")
    plt.title("Hist angular difference between true orientation vector and fit orientation\n on Random Orientations Data")
    plt.tight_layout()
    plt.savefig("grid_random_orientation_errors.png")
    plt.clf()
    print("Matched Position of Every Test Particle on lattice with random"
          " orientations\n Max angle between true orientation vector and fit"
          " orientation vector: ", np.max(angle_diff), " degrees")



def test_random(tv_fname='grid_random_la.npy',
                fit_name='grid_random/grid_random0.5_positions_orientation.npy'
                ):
    true_values = np.load(tv_fname, allow_pickle=True)
    fit_values = np.load(fit_name, allow_pickle=True)
    
    fig, ax = plt.subplots(2)
    position_diff = []
    angle_diff = []
    
    for i in true_values:
        distance = 1000
        match = [None]
        for j in fit_values:
            x = np.sum(np.abs(i[0] - j[0]))
            
            if x < distance:
                distance = x
                match = j
        
        position_diff.append(distance)
        
        if distance > 10:
            print("Likely did not match all particles", distance)
    
            print(i[0])
            print()
            print(match[0])
        else:
            angle_diff.append(
                np.rad2deg(np.arccos(np.abs(np.dot(i[2], match[2])))))

    ax[0].hist(position_diff, bins=30)
    ax[0].set_xlabel("Voxels")
    ax[0].set_ylabel("Count")
    ax[0].set_title(
        "Hist center of mass difference between true orientation "
        "and fit orientation\n"
        " on Random Orientations and position Data")
    
    ax[1].hist(angle_diff, bins=30)
    ax[1].set_xlabel("Degrees")
    ax[1].set_ylabel("Count")
    ax[1].set_title("Hist angular difference between true "
                    "orientation and fit orientation\n"
                    "on Random Orientations and position Data")
    plt.tight_layout()
    plt.savefig("grid_random_errors.png")
    plt.clf()

    print("Matched Position of Every Test Particle on with random"
          " position and orientations\n Max angle between true "
          "orientation vector and fit"
          " orientation vector: ", np.max(angle_diff), " degrees\n"
          "Max distance between fit center of mass and true center of mass ",
          np.max(position_diff), " voxels\n"
          "Mean Errors with stdev:\n"
          "Angle: ", np.mean(angle_diff), " +/- ", np.sqrt(np.var(angle_diff)),
          "\nPosition: ", np.mean(position_diff), " +/- ",
          np.sqrt(np.var(position_diff)),
          "\nNote this is for BINARIZED particles with characteristic scale 25"
          "\nvoxels and will vary for smaller and larger ones, as well as with "
          "\nnoise and artifacts from the tomography scan process and how much "
          "\nerosion takes place at the binirize step."

          )


def main():
    test_grid()
    test_random_orientations()
    test_random()


if __name__ == '__main__':
    main()

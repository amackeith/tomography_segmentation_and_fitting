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
    true_values = np.load('../old_py/grid_lattice_la.npy', allow_pickle=True)
    fit_values = np.load('../old_py/grid_lattice/grid_lattice_positions_orientation.npy',
                         allow_pickle=True)
    
    angle_diff = []
    for i in true_values:
        distance = 1000
        match = None
        for j in fit_values:
            x = np.sum(np.abs(i[0]+200 - j[0]))
            
            if x < distance:
                distance = x
                match = j
        
        if distance != 0:
            print("Failed to match every test particle", distance)
            print(i)
            print(match)
            exit()
            
        else:
            angle_diff.append(
                np.rad2deg(np.arccos(np.abs(np.dot(i[2], match[2])))))
    
    print("Matched Position of Every Test Particle")


def test_random_orientations():
    true_values = np.load('../old_py/grid_lattice_random_orientations_la.npy', allow_pickle=True)
    fit_values = np.load(
        '../old_py/grid_lattice_random_orientations/grid_lattice_random_orientations_positions_orientation.npy',
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
            exit()
        else:
            angle_diff.append(
                np.rad2deg(np.arccos(np.abs(np.dot(i[2], match[2])))))
    
    plt.hist(angle_diff, bins=30)
    plt.xlabel("Degrees")
    plt.ylabel("Count")
    plt.title("Hist angular difference between true orientation vector and fit orientation\n on Random Orientations Data")
    plt.tight_layout()
    plt.show()


def test_random():
    true_values = np.load('../old_py/grid_random_la.npy', allow_pickle=True)
    fit_values = np.load('../old_py/grid_random/grid_random_positions_orientation.npy',
                         allow_pickle=True)
    
    fig, ax = plt.subplots(2)
    position_diff = []
    angle_diff = []
    
    for i in true_values:
        distance = 1000
        match = None
        for j in fit_values:
            x = np.sum(np.abs(i[0] - j[0]))
            
            if x < distance:
                distance = x
                match = j
        
        position_diff.append(distance)
        
        if distance > 0.2:
            print(i[0])
            print()
            print(match[0])
            print("Failed to match every test particle", distance)
            exit()
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
    plt.show()


def main():
    test_grid()
    test_random_orientations()
    test_random()


if __name__ == '__main__':
    main()

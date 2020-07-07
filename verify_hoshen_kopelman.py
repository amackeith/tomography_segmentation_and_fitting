import numpy as np
import matplotlib.pyplot as plt

from watershedtools_py import \
    hoshen_kopelmann3d, better_labels, hoshen_kopelmann_verify_3d

# I am leaving this test code so that others may see how the code works
# on a simple problem
# take a look at arr, then also hoshen_kopelmann3d(arr) to get an idea.
arr = np.array([[1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                [0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1],
                [0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1],
                [1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
                [0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1],
                [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
                [1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                [1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1],
                [1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0]])

arr2 = np.array([[1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                 [0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1],
                 [0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                 [1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
                 [0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1],
                 [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
                 [1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                 [1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1],
                 [1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
                 [0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0]])
arr0 = np.zeros(arr.shape)
arr = np.array(
    [arr, arr, arr2, arr2, arr2, arr2, arr0, arr2,
     arr2, arr2, arr2, arr2, arr2, arr, arr, arr2, arr2, arr2, arr2, arr0,
     arr2, arr2, arr2, arr2, arr2, arr2])

arr = np.array([arr[20], arr[20]])

arr = np.load("hoshen_testing.npy")
arr = hoshen_kopelmann3d(arr)
arr = better_labels(arr)
# print(set(labels))
# plt.imshow(arr[0])
# plt.show()
hoshen_kopelmann_verify_3d(arr)

for i in range(arr.shape[0]):
    if i % 4 != 0:
        continue
    arr2 = arr.copy()
    arr2 = arr2.flatten()
    arr2 = list(set(arr2))
    print(len(arr2))
    plt.imshow(arr[i])
    plt.title(str(i))
    plt.show()

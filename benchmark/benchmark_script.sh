#!/bin/bash


python create_benchmark_volumes.py

tomography_segment_and_fit.py -i grid_lattice.npy -o grid_lattice -de_noised_volume -b 0.5

tomography_segment_and_fit.py -i grid_lattice_random_orientations.npy -o grid_lattice_random_orientations -de_noised_volume  -b 0.5

tomography_segment_and_fit.py -i grid_random.npy -o grid_random -de_noised_volume  -b 0.5

python compare_benchmark.py


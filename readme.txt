This software is for segmentation and fitting of 3D density data from
tomographic scans of pickings of spheriods.

Use:
This takes in a tomography as input and outputs a list of

Start by installing python dependencies.

Since tomopy is only availible through conda-forge the most straight forward way to do
this (maybe the only way) is through conda. To do this install conda and then
in the directory you downloaded the repo to run:

$ conda env create -f seg_fit_conda_environment.yml

then activate the environment with

$ conda activate seg_fit

to add tomography_segment_and_fit.py to the path:
on mac:

$ echo "export PATH=\$PATH:$(pwd)" >> ~/.bash_profile

or on linux:

$ echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc


Install swig:
linux:
sudo apt-get install swig

In order to get the cpp on Linux or mac:
$ cd watershedtools_cpp
$ bash ./build.sh

the cpp should work on windows as well, 
but you will need to figure out how to build it
yourself.

now you should be able to run this with:

$ tomography_segment_and_fit.py -i input_file.npy -o output_folder -other_options

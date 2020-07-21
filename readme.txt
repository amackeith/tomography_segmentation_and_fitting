About:
This software was written to segment (and fit for position and orientation)
 tomography scans of granular samples of spheroids.


Start by installing python dependencies.

Since tomopy is only availible through conda-forge the most straight forward way to do
this (maybe the only way) is through conda. To do this install conda and then
in the directory you downloaded the repo to run:

$ conda env create -f seg_fit_conda_environment.yml

then activate the environment with (note this env must be active when you call
tomography_segment_and_fit.py from the command line in other directories)

$ conda activate seg_fit


to add tomography_segment_and_fit.py to the path:
On OS X:

$ echo "export PATH=\$PATH:$(pwd)" >> ~/.bash_profile

or on linux:

$ echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc


##################################
CPP acceleration
##################################
(The code can be run without this using the python version of functions
but it is significantly slower).


Install swig:
linux:
sudo apt-get install swig

Or for other platforms (though mac may come with it):
https://www.dev2qa.com/how-to-install-swig-on-macos-linux-and-windows/

In order to get the cpp on Linux or mac:
$ cd watershedtools_cpp
$ bash ./build.sh

the cpp should work on windows as well, 
but you will need to figure out how to build it
yourself.

now you should be able to run this with (where input_file.npy) is in that directory.

$ tomography_segment_and_fit.py -i input_file.npy -o /path/to/output_folder -other_options

for more info run
$ tomography_segment_and_fit.py -h


Usage: to test if setup has worked run:
$ cd benchmark
$ ./benchmark_script.sh

benchmark_script.sh provides simple example usage.
It creates three digitally fabricated "packings" and then checks the accuracy
of the fit on them. You should get results similar to those displayed in the folder
or exactly the same depending on the seed.

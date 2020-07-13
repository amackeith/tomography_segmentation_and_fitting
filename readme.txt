Start by installing python dependencies.

Since tomopy is only availible through conda-forge the most straight forward way to do
this (maybe the only way) is through conda. To do this install conda and run:

conda env create -f seg_fit_conda environment.yml

then activate the environment with
conda activate seg_fit 


In order to get the cpp on Linux or mac:
cd watershedtools_cpp
bash ./build.sh

the cpp should work on windows as well, 
but you will need to figure out how to build it
yourself.


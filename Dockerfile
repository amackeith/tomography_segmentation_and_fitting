FROM ubuntu:18.04




RUN apt update && apt install -y swig
RUN apt install -y wget vim tree 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda
ENV PATH=$PATH:/miniconda/condabin:/miniconda/bin

# update and install dependencies
RUN         apt update && \
            apt install -y build-essential 
#libpython3.6-dev python3.6



WORKDIR /home/workdir
COPY seg_fit_conda_environment.yml /home/workdir


## set up conda env
RUN conda env create -f seg_fit_conda_environment.yml
RUN conda init bash
RUN echo "conda activate seg_fit" >> ~/.bashrc
ENV PATH="${PATH}:/home/workdir"
COPY . /home/workdir/

####### build C++ code

# compile the C++ source code  -fPIC flag is needed for what we want to do
RUN cd watershedtools_cpp && g++ --std=c++11 -fPIC -c watershedtools_cpp.cpp

# use SWIG on the interface file 'watershedtools_cpp.i'
# this creates the wrapper
# code ('watershedtools_cpp_wrap.cpp'):

RUN  cd watershedtools_cpp && swig -python -c++ -o watershedtools_cpp_wrap.cpp watershedtools_cpp.i

# then compile the wrapper, this step must link to the Python.h header as
# well as the numpy header.
# NOTE: you must have set up and activated the conda env before this step.
RUN  cd watershedtools_cpp && g++ -std=c++11 -fPIC -c \
      -I /miniconda/envs/seg_fit/include/python3.6m \
      -I /miniconda/envs/seg_fit/lib/python3.6/site-packages/numpy/core/include \
        watershedtools_cpp_wrap.cpp
# Link both object files into a shared library '_watershedtools_cpp.so'.
# Python will search for it under this name: watershedtools_cpp

RUN  cd watershedtools_cpp &&  g++ --std=c++11 -shared watershedtools_cpp.o watershedtools_cpp_wrap.o -o \
      _watershedtools_cpp.so

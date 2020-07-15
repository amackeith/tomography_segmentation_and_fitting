#!/bin/bash

# Build Python module watershedtools_cpp:
# NOTE
# this greates wrapper code with c++, but gcc could work as well and maybe ld


# compile the C++ source code  -fPIC flag is needed for what we want to do

gcc --std=c++11 -fPIC -c watershedtools_cpp.cpp

# use SWIG on the interface file 'watershedtools_cpp.i'
# this creates the wrapper
# code ('watershedtools_cpp_wrap.cpp'):

swig -python -c++ -o watershedtools_cpp_wrap.cpp watershedtools_cpp.i

# then compile the wrapper, this step must link to the Python.h header as
# well as the numpy header.
# NOTE: you must have set up and activated the conda env before this step.
gcc -std=c++11 -fPIC -c \
      -I "${CONDA_PREFIX}"/include/python3.6m \
      -I "${CONDA_PREFIX}"/lib/python3.6/site-packages/numpy/core/include \
        watershedtools_cpp_wrap.cpp
# Link both object files into a shared library '_watershedtools_cpp.so'.
# Python will search for it under this name: watershedtools_cpp

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
gcc --std=c++11 -shared watershedtools_cpp.o watershedtools_cpp_wrap.o -o \
      _watershedtools_cpp.so -lm
elif [[ "$OSTYPE" == "darwin"* ]]; then
ld -bundle -flat_namespace -undefined suppress -o _watershedtools_cpp.so \
      watershedtools_cpp.o watershedtools_cpp_wrap.o
else
  raise error "bash script compile of CPP only supported on Linux and Mac OS, the python version will still work though"
fi

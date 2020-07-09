#!/bin/bash

# Build Python module watershedtools_cpp:
# NOTE
# this greates wrapper code with c++, but gcc could work as well and maybe ld


# compile the C++ source code  -fPIC flag is needed for what we want to do

g++ -fPIC -c watershedtools_cpp.cpp

# use SWIG on the interface file 'watershedtools_cpp.i'
# this creates the wrapper
# code ('watershedtools_cpp_wrap.cpp'):

swig -python -c++ -o watershedtools_cpp_wrap.cpp watershedtools_cpp.i

# then compile the wrapper
# depending on where your Python.h file is (and your version of python).

g++ -fPIC -c -I /usr/include/python3.6 watershedtools_cpp_wrap.cpp

# Link both object files into a shared library '_watershedtools_cpp.so'.
# Python will search for it under this name: watershedtools_cpp

g++ -shared watershedtools_cpp.o watershedtools_cpp_wrap.o -o \
      _watershedtools_cpp.so -lm

In order to use the cpp acceleration for this project
run ./build.sh

you may need to use a different version of numpy.i that is compatible with your
version of python / numpy. 

You may also need to change the line 
g++ -fPIC -c -I /usr/include/python3.6 watershedtools_cpp_wrap.cpp

to reflect where your copy of Python.h is.  If you do not have it (on linux)
you need to install python.X.x-dev where X.x is you version of python e.g 3.8
or 3.6. (The rest of the code is not compatible with 2.7 so not that). 

so: 
sudo apt-get install python3.6-dev or something along those lines.



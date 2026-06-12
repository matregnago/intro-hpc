#!/bin/sh
# whatis: standard development environment to use the chameleon library

sudo apt-get update -y
sudo apt-get install build-essential git gfortran cmake doxygen python-is-python3 pkg-config libopenblas-dev liblapacke-dev libstarpu-dev libopenmpi-dev libhwloc-dev -y
# libopenblas-dev can be replaced by libmkl-dev or liblapack-dev

# install chameleon with cmake
git clone --recursive https://gitlab.inria.fr/solverstack/chameleon.git
cd chameleon && cmake -S . -B build -DBUILD_SHARED_LIBS=ON && cmake --build build -j5 && sudo cmake --install build
# enable MPI with -DCHAMELEON_USE_MPI=ON

# example usage: use chameleon drivers for performance tests
chameleon_stesting -t 4 -o gemm -n 2000
# if MPI enabled
mpiexec -n 2 chameleon_stesting -t 4 -o gemm -n 4000

# example usage: use chameleon library in your own cmake project (we provide a CHAMELEONConfig.cmake)
cd ..
git clone https://gitlab.inria.fr/solverstack/distrib.git
cd distrib/cmake/test/chameleon && mkdir build && cd build && cmake .. && make && ./test_chameleon

# example usage: use chameleon library in your own not cmake project
# use pkg-config to get compiler flags and linking
pkg-config --cflags chameleon
pkg-config --libs chameleon
# if there are static libraries use the --static option of pkg-config

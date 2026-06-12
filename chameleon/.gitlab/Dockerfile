#
#  @file Dockerfile
#
#  @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @date 2025-12-19
#
# This docker image is used to test the package during gitlab-ci pipelines.
# It should be stored in the gitlab's project container registry:
# https://gitlab.inria.fr/solverstack/chameleon/container_registry
#
FROM registry.gitlab.inria.fr/solverstack/docker/analysis

USER root

# Chameleon dependencies: cpu kernels
RUN apt-get -y update && \
    apt-get -y upgrade --no-install-recommends \
    libblas-dev \
    liblapack-dev \
    liblapacke-dev
RUN update-alternatives --set libblas.so-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/blas/libblas.so
RUN update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/blas/libblas.so.3
RUN update-alternatives --set liblapack.so-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/lapack/liblapack.so
RUN update-alternatives --set liblapack.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3

# Starpu dependencies
RUN apt-get -y upgrade --no-install-recommends \
    autoconf \
    automake \
    libtool-bin \
    libfxt-dev \
    libhwloc-dev \
    libopenmpi-dev \
    libsimgrid-dev

# To enable eztrace chameleon module
RUN apt-get -y upgrade --no-install-recommends \
    eztrace \
    libeztrace-dev \
    libotf2-trace-dev

RUN apt-get autoremove -y
RUN apt-get autoclean -y
RUN apt-get purge -y

RUN update-ca-certificates

# Fix eztrace.pc
RUN sed -i 's+/usr//usr/+/usr/+' /usr/lib/x86_64-linux-gnu/pkgconfig/eztrace.pc

# Install PARSEC
RUN git clone https://bitbucket.org/mfaverge/parsec.git && \
    cd parsec && \
    git checkout mymaster && \
    git checkout 6022a61dc96c25f11dd2aeabff2a5b3d7bce867d && \
    git submodule update && \
    cmake -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_DIST_WITH_MPI=OFF && \
    cmake --build build -j5 && \
    cmake --install build && \
    cd .. && \
    rm parsec/ -rf

# Install QUARK
RUN git clone https://github.com/ecrc/quark && \
    cd quark/ && \
    git checkout db4aef9a66a00487d849cf8591927dcebe18ef2f && \
    sed -i -e "s#prefix=.*#prefix=/usr#g" make.inc && \
    sed -i -e "s#CFLAGS=.*#CFLAGS= -g -O2 -DADD_ -fPIC#g" make.inc && \
    make && \
    make install && \
    cd .. && \
    rm quark/ -rf

# Install StarPU
RUN wget https://files.inria.fr/starpu/starpu-1.4.12/starpu-1.4.12.tar.gz > /dev/null 2>&1 && \
    tar xzf starpu-1.4.12.tar.gz && \
    cd starpu-1.4.12/ && \
    ./configure --prefix=/home/gitlab/install/starpu --enable-debug --disable-build-doc --with-fxt --disable-simgrid && \
    make -j5 install

# Install StarPU with Simgrid
RUN cd starpu-1.4.12/ && \
    ./configure --prefix=/home/gitlab/install/starpu-simgrid --enable-debug --disable-build-doc --with-fxt --enable-simgrid && \
    make -j5 install && \
    cd ../ && \
    rm starpu-1.4.12/ -rf && \
    rm starpu-1.4.12.tar.gz

USER gitlab

# Set some default environment variables
ENV LIBRARY_PATH=/usr/local/lib
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV PKG_CONFIG_PATH=/home/gitlab/install/starpu/lib/pkgconfig

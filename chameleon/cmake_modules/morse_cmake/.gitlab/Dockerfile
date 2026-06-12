FROM ubuntu:24.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update

RUN apt-get -y upgrade --no-install-recommends \
    # linux base packages
    bash ca-certificates cpio sudo vim wget \
    # dev tools
    bison cppcheck flex git pkg-config python-is-python3 \
    # C/Fortran compiling tools
    build-essential cmake cmake-data cmake-curses-gui gfortran

RUN apt-get -y upgrade --no-install-recommends \
    # packages to find with morse_cmake
    libblas-dev \
    libfftw3-dev \
    libfxt-dev \
    liblapack-dev \
    liblapacke-dev \
    libmetis-dev \
    libmumps-dev \
    libopenmpi-dev \
    libpapi-dev \
    libparmetis-dev \
    libptscotch-dev \
    libscalapack-mpi-dev \
    libsimgrid-dev \
    libstarpu-dev \
    libsuitesparse-dev \
    petsc-dev \
    slepc-dev

RUN apt-get autoremove -y
RUN apt-get autoclean -y
RUN apt-get purge -y

RUN update-ca-certificates
RUN update-alternatives --set libblas.so-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/blas/libblas.so
RUN update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/blas/libblas.so.3
RUN update-alternatives --set liblapack.so-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/lapack/liblapack.so
RUN update-alternatives --set liblapack.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3

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

# Install PARSEC
RUN git clone https://bitbucket.org/mfaverge/parsec.git && \
    cd parsec && \
    git checkout mymaster && \
    git checkout 6022a61dc96c25f11dd2aeabff2a5b3d7bce867d && \
    git submodule update && \
    cmake -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_DIST_WITH_MPI=ON && \
    cmake --build build -j5 && \
    cmake --install build && \
    cd .. && \
    rm parsec/ -rf

RUN groupadd -f -g 1001 gitlab && \
    useradd -u 1001 -g gitlab -ms /bin/bash gitlab && \
    echo "gitlab:gitlab" | chpasswd && adduser gitlab sudo

USER gitlab

# change the default shell to be bash
SHELL ["/bin/bash", "-c"]

# default working directory is
WORKDIR /home/gitlab

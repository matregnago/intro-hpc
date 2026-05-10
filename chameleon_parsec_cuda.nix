{
  stdenv,
  lib,
  fetchgit,
  cmake,
  pkg-config,
  gfortran,
  openblas,
  lapack,
  hwloc,
  openmpi,
  parsec,
  python3,
  cudaPackages,
}:

stdenv.mkDerivation {
  pname = "chameleon-parsec-cuda";
  version = "1.4.0";

  src = fetchgit {
    url = "https://gitlab.inria.fr/solverstack/chameleon/";
    rev = "a907fde133b122b51418f04528d15fa387b6ad28";
    fetchSubmodules = true;
    hash = "sha256-qqKrdzybcFVtmcp+2gVECChvEuD53rVSBT1LTflRuQY=";
  };

  nativeBuildInputs = [
    cmake
    pkg-config
    gfortran
    python3
    # nvcc must be in PATH so CMake's check_language(CUDA) finds it
    cudaPackages.cuda_nvcc
  ];

  buildInputs = [
    openblas
    lapack
    hwloc
    openmpi
    parsec
    cudaPackages.cuda_cudart  # CUDA::cudart
    cudaPackages.cuda_cccl    # nv/target (included by cuda_fp16.h in CUDA 12+)
    cudaPackages.libcublas    # CUDA::cublas
    cudaPackages.libcusolver  # CUDA::cusolver
  ];

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
    "-DCHAMELEON_SCHED=PARSEC"
    "-DCHAMELEON_USE_MPI=OFF"
    "-DCHAMELEON_USE_CUDA=ON"
    "-DCMAKE_CUDA_COMPILER=${lib.getExe' cudaPackages.cuda_nvcc "nvcc"}"
  ];
}

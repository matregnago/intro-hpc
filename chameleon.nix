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
  starpu,
  parsec,
  python3,
  cudaPackages ? null,
  autoAddDriverRunpath ? null,
  enableCuda ? true,
  runtime ? "starpu",
}:

let
  schedMap = {
    starpu = "STARPU";
    parsec = "PARSEC";
  };
  sched = schedMap.${runtime};
in
stdenv.mkDerivation {
  pname = "chameleon-${runtime}";
  version = "1.4.0";

  src = fetchgit {
    url = "https://gitlab.inria.fr/solverstack/chameleon.git";
    rev = "13960d0d8e529777bd7242d9f4563253b81dcb3b";
    fetchSubmodules = true;
    hash = "sha256-vUTdzbaLn0eYgQU7TDdSRknFFV0MEI6h6o8DYJfP7y4=";
  };

  patches = [ ./patches/chameleon.patch ];

  nativeBuildInputs = [
    cmake
    pkg-config
    gfortran
    python3
  ]
  ++ lib.optionals enableCuda [
    cudaPackages.cuda_nvcc
    autoAddDriverRunpath
  ];

  buildInputs = [
    openblas
    lapack
    hwloc
    openmpi
  ]
  ++ lib.optionals enableCuda [
    cudaPackages.cuda_cudart # CUDA::cudart
    cudaPackages.cuda_cccl # nv/target (included by cuda_fp16.h in CUDA 12+)
    cudaPackages.cuda_nvml_dev # nvml.h (starpu_cuda.h includes it when StarPU was built with NVML)
    cudaPackages.libcublas # CUDA::cublas
    cudaPackages.libcusolver # CUDA::cusolver
    cudaPackages.libcusparse # cusparse.h (starpu_cusparse.h includes it when StarPU was built with cuSPARSE)
  ]
  ++ lib.optional (runtime == "starpu") starpu
  ++ lib.optional (runtime == "parsec") parsec;

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
    "-DCHAMELEON_SCHED=${sched}"
    "-DCHAMELEON_USE_MPI=OFF"
  ]
  ++ (
    if enableCuda then
      [
        "-DCHAMELEON_USE_CUDA=ON"
        "-DCMAKE_CUDA_COMPILER=${lib.getExe' cudaPackages.cuda_nvcc "nvcc"}"
      ]
    else
      [
        "-DCHAMELEON_USE_CUDA=OFF"
      ]
  );
}

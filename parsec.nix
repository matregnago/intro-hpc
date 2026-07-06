{
  stdenv,
  lib,
  fetchgit,
  cmake,
  pkg-config,
  hwloc,
  gfortran,
  flex,
  bison,
  python3,
  fetchFromGitHub,
  autoconf,
  automake,
  libtool,
  zlib,
  cudaPackages,
  autoAddDriverRunpath,
}:

let
  gtg = stdenv.mkDerivation {
    pname = "gtg";
    version = "unstable";
    src = fetchFromGitHub {
      owner = "trahay";
      repo = "Generic-Trace-Generator";
      rev = "master";
      sha256 = "sha256-lc4tPaPE7Q4HnLlH2z3FaxXCX7oXFIUJszcpKlAeOFg=";
    };
    nativeBuildInputs = [
      autoconf
      automake
      libtool
      pkg-config
      gfortran
    ];
    buildInputs = [ zlib ];
    preConfigure = ''
      ./autogen.sh || (autoreconf -fi)
    '';
  };
in
stdenv.mkDerivation {
  pname = "parsec";
  version = "mymaster";

  src = fetchgit {
    url = "https://bitbucket.org/mfaverge/parsec.git";
    rev = "83ff82ec9176513c912d46114ca6ae211a2720ee";
    hash = "sha256-Esp8jUM3OLtliLbsnIR1jCViEbPzZANDt9fsZkj2hys=";
  };

  nativeBuildInputs = [
    cmake
    pkg-config
    gfortran
    flex
    bison
    python3
    cudaPackages.cuda_nvcc
    autoAddDriverRunpath
  ];

  buildInputs = [
    hwloc
    stdenv.cc.cc.lib
    gtg
    cudaPackages.cuda_cudart
  ];
  patches = [ ./patches/parsec.patch ];

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
    "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
    "-DPARSEC_DIST_WITH_MPI=OFF"
    "-DPARSEC_PROF_TRACE=ON"
    "-DPINS_ENABLE=ON"
    "-DPARSEC_PROF_GRAPHER=ON"
    "-DPARSEC_PROF_TRACE_SCHEDULING_EVENTS=ON"
    "-DPARSEC_PROF_TRACE_ACTIVE_ARENA_SET=ON"
    "-DPARSEC_GPU_WITH_CUDA=ON"
    "-DCMAKE_CUDA_COMPILER=${lib.getExe' cudaPackages.cuda_nvcc "nvcc"}"
  ];
  env.NIX_CFLAGS_COMPILE = toString [
    "-Wno-error=implicit-function-declaration"
    "-Wno-error=incompatible-pointer-types"
  ];
}

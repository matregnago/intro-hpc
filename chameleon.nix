{
  stdenv,
  fetchurl,
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
  eztrace,
}:

stdenv.mkDerivation rec {
  pname = "chameleon";
  version = "1.4.0";

  src = fetchurl {
    url = "https://gitlab.inria.fr/api/v4/projects/616/packages/generic/source/v${version}/chameleon-${version}.tar.gz";
    hash = "sha256-e2WmFo69/BGhyZlEVO5HRnP8kjSNWqhyToPDT9UBCq8=";
  };

  nativeBuildInputs = [
    cmake
    pkg-config
    gfortran
    python3
  ];

  buildInputs = [
    openblas
    lapack
    hwloc
    openmpi
    starpu
    parsec
    eztrace
  ];

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
    "-DCHAMELEON_SCHED=PARSEC"
    "-DCHAMELEON_USE_MPI=OFF"
    "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
    "-DCHAMELEON_KERNELS_TRACE=ON"
  ];
}

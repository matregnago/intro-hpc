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
  eztrace,
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
    url = "https://gitlab.inria.fr/solverstack/chameleon/";
    rev = "a907fde133b122b51418f04528d15fa387b6ad28";
    fetchSubmodules = true;
    hash = "sha256-qqKrdzybcFVtmcp+2gVECChvEuD53rVSBT1LTflRuQY=";
  };

  # patches = [ ./chameleon.patch ];
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
    eztrace
  ]
  ++ lib.optional (runtime == "starpu") starpu
  ++ lib.optional (runtime == "parsec") parsec;

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
    "-DCHAMELEON_SCHED=${sched}"
    "-DCHAMELEON_USE_MPI=OFF"
  ];
}

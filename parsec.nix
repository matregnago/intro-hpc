{
  stdenv,
  fetchgit,
  cmake,
  pkg-config,
  hwloc,
  openmpi,
  gfortran,
  flex,
  bison,
  python3,
  fetchFromGitHub,
  autoconf,
  automake,
  libtool,
  zlib,
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
    rev = "6022a61dc96c25f11dd2aeabff2a5b3d7bce867d";
    fetchSubmodules = true;
    hash = "sha256-iJEIJLdakIxGFz+qXNTcbGOISzXbI/XLVq45k4GqhoM=";
  };

  nativeBuildInputs = [
    cmake
    pkg-config
    gfortran
    flex
    bison
    python3
  ];

  buildInputs = [
    hwloc
    openmpi
    stdenv.cc.cc.lib
    gtg
  ];
  patches = [ ./parsec.patch ];

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
    "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
    "-DPARSEC_PROF_TRACE=ON"
    "-DPINS_ENABLE=ON"
    "-DPARSEC_PROF_GRAPHER=ON"
    "-DPARSEC_PROF_TRACE_SCHEDULING_EVENTS=ON"
    "-DPARSEC_PROF_TRACE_ACTIVE_ARENA_SET=ON"
  ];
  env.NIX_CFLAGS_COMPILE = toString [
    "-Wno-error=implicit-function-declaration"
    "-Wno-error=incompatible-pointer-types"
  ];
}
